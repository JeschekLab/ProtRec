### Setup
# Import Python modules
# Plotting is done with plotnine, which uses ggplot2-like syntax.
# The analysis uses the Bloom lab's alignparse and dms_variants packages.
import collections
import math
import os
import re
import time
import warnings

import alignparse
import alignparse.ccs
from alignparse.constants import CBPALETTE
import alignparse.minimap2
import alignparse.targets

import dms_variants
import dms_variants.plotnine_themes
import dms_variants.utils

from IPython.display import display, HTML

import numpy

import pandas as pd

from plotnine import *

import yaml

# Set plotnine theme to the one defined in dms_variants:
theme_set(dms_variants.plotnine_themes.theme_graygrid())

# Versions of key software:
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")

# Ignore warnings that clutter output:
warnings.simplefilter('ignore')

# Read the configuration file:
with open('config.yaml') as f:
    config = yaml.safe_load(f)

# Make output directory for figures:
os.makedirs(config['figs_dir'], exist_ok=True)
os.makedirs(config['process_ccs_dir'], exist_ok=True)

### PacBio amplicons

# Get the amplicons sequenced by PacBio as the alignment target along with the
# specs on how to parse the features:
print(f"Reading amplicons from {config['amplicons']}")
print(f"Reading feature parse specs from {config['feature_parse_specs']}")

targets = alignparse.targets.Targets(
                seqsfile=config['amplicons'],
                feature_parse_specs=config['feature_parse_specs'],
                allow_clipped_muts_seqs=True)

# Draw the target amplicons:
fig = targets.plot(ax_width=7,
                   plots_indexing='biopython',  # numbering starts at 0
                   ax_height=2,  # height of each plot
                   hspace=1.2,  # vertical space between plots
                   )

plotfile = os.path.join(config['figs_dir'], 'amplicons.pdf')
print(f"Saving plot to {plotfile}")
fig.savefig(plotfile, bbox_inches='tight')

# Write out the specs used to parse the features (these are the same specs
# provided as feature_parse_specs when initializing targets, but with defaults
# filled in):
print(targets.feature_parse_specs('yaml'))

### CCS stats for PacBio runs

# Set data frame with information on PacBio runs:
pacbio_run = pd.DataFrame(
    data = {
    'name': [config['project_name']],
    'fastq': [config['pacbio_run']]})

# Create an object that summarizes the ccs runs:
ccs_summaries = alignparse.ccs.Summaries(pacbio_run,
                                         report_col=None,
                                         ncpus=config['max_cpus'],
                                         )

# If available, plot statistics on the number of ZMWs for each run:
if ccs_summaries.has_zmw_stats():
    p = ccs_summaries.plot_zmw_stats()
    p = p + theme(panel_grid_major_x=element_blank()) # no vert. grid lines
    _ = p.draw()
else:
    print('No ZMW stats available.')

# # Plot statistics on generated CCSs: their length, number of subread passes,
# # and accuracy (as reported by the ccs program):
# for variable in ['length', 'passes', 'accuracy']:
#     if ccs_summaries.has_stat(variable):
#         p = ccs_summaries.plot_ccs_stats(variable, maxcol=7, bins=25)
#         p = p + theme(panel_grid_major_x=element_blank()) # no vert. grid lines
#         _ = p.draw()
#     else:
#         print(f"No {variable} statistics available.")


### Align CCSs to amplicons

# We now align the CCSs to the amplicon and parse features from the resulting
# alignments using the specs above.
# First, we initialize an alignparse.minimap2.Mapper to align the reads to SAM
# files:
mapper = alignparse.minimap2.Mapper(alignparse.minimap2.OPTIONS_CODON_DMS)

print(f"Using `minimap2` {mapper.version} with these options:\n" +
      ' '.join(mapper.options))

# Next, we use Targets.align_and_parse to create the alignments and parse them:
readstats, aligned, filtered = targets.align_and_parse(
        df=pacbio_run,
        mapper=mapper,
        outdir=config['process_ccs_dir'],
        # name_col='run',
        # group_cols=['name', 'library'],
        queryfile_col='fastq',
        overwrite=True,
        ncpus=config['max_cpus'],
        )

# First, examine the read stats from the alignment / parsing, both extracting
# alignment target name and getting stats aggregated by target:
readstats = (
    readstats
    .assign(category_all_targets=lambda x: x['category'].str.split().str[0],
            target=lambda x: x['category'].str.split(None, 1).str[1],
            valid=lambda x: x['category_all_targets'] == 'aligned')
    )

# Plot the read stats:
ncol = 7
# p = (
#     ggplot(readstats
#            .groupby(['name', 'category_all_targets', 'valid'])
#            .aggregate({'count': 'sum'})
#            .reset_index(),
#            aes('category_all_targets', 'count', fill='valid')) +
#     geom_bar(stat='identity') +
#     facet_wrap('~ name', ncol=ncol) +
#     theme(axis_text_x=element_text(angle=90),
#           figure_size=(1.85 * min(ncol, len(pacbio_run)),
#                        2 * math.ceil(len(pacbio_run) / ncol)),
#           panel_grid_major_x=element_blank(),
#           legend_position='none',
#           ) +
#     scale_fill_manual(values=CBPALETTE)
#     )
# _ = p.draw()

# Now let's see why we filtered the reads. First, we do some transformations on
# the filtered dict returned by Targets.align_and_parse. Then we count up the
# number of CCSs filtered for each reason, and group together "unusual" reasons
# that represent less than some fraction of all filtering. For now, we group
# together all targets to the stats represent all targets combined:
other_cutoff = 0.02  # group as "other" reasons with <= this frac

filtered_df = (
    pd.concat(df.assign(target=target) for target, df in filtered.items())
    .groupby(['name', 'filter_reason'])
    .size()
    .rename('count')
    .reset_index()
    .assign(tot_reason_frac=lambda x: (x.groupby('filter_reason')['count']
                                       .transform('sum')) / x['count'].sum(),
            filter_reason=lambda x: numpy.where(x['tot_reason_frac'] > other_cutoff,
                                                x['filter_reason'],
                                                'other')
            )
    )


# Now plot the filtering reason for all runs:
ncol = 7
nreasons = filtered_df['filter_reason'].nunique()

# p = (
#     ggplot(filtered_df, aes('filter_reason', 'count')) +
#     geom_bar(stat='identity') +
#     facet_wrap('~ name', ncol=ncol) +
#     theme(axis_text_x=element_text(angle=90),
#           figure_size=(0.25 * nreasons * min(ncol, len(pacbio_run)),
#                        2 * math.ceil(len(pacbio_run) / ncol)),
#           panel_grid_major_x=element_blank(),
#           )
#     )
# _ = p.draw()

# Finally, we take the successfully parsed alignments and read them into a data
# frame, keeping track of the target that each CCS aligns to. We also drop the
# pieces of information we won't use going forward, and rename a few columns:
aligned_df = (
    pd.concat(df.assign(target=target) for target, df in aligned.items())
    .drop(columns=['query_clip5', 'query_clip3', 'name'])
    .rename(columns={'barcode_sequence': 'barcode'})
    )

### Write valid CCSs

# Write the processed CCSs to a file:

aligned_df.to_csv(config['processed_ccs_file'], index=False)

print("Barcodes and mutations for valid processed CCSs "
      f"have been written to {config['processed_ccs_file']}.")

# In the next notebook, we analyze these processed CCSs to build the variants.
