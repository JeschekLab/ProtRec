# tested
# 3.9.15 not working


# install stuff
conda install python=3.9.13
conda install -c bioconda --force-reinstall minimap2

pip install --force-reinstall alignparse
pip install --force-reinstall dms_variants
# pip install --force-reinstall dna_features_viewer

# utils
conda install -c conda-forge --force-reinstall urllib3
conda install -c anaconda --force-reinstall charset-normalizer
conda update --all
