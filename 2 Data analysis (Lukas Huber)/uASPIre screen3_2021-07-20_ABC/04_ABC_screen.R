# load libraries
library(tidyverse)
library(viridis)
library(cowplot)

#define which package should be preferred for certain functions
library(conflicted)
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")

# create folders ----------------------------------------------------------

# create a folder to store the individual plots and data in
dir.create("./04_ABC_screen")
dir.create("./04_ABC_screen/results")
dir.create("./04_ABC_screen/plots")
dir.create("./04_ABC_screen/plots/heatmaps")
dir.create("./04_ABC_screen/plots/heatmaps/orthogonal")
dir.create("./04_ABC_screen/plots/heatmaps/promiscuous")

# prepare data ------------------------------------------------------------

# read data
rep1 <- read.table("./03_false_positives_replicate_1/results/Data_ABC_fc_corrected_min20proteinlevel.txt", sep = "\t", header = TRUE) %>%
  mutate(rep = 1)
rep2 <- read.table("./03_false_positives_replicate_2/results/Data_ABC_fc_corrected_min20proteinlevel.txt", sep = "\t", header = TRUE) %>%
  mutate(rep = 2)

data <- rbind(rep1, rep2)

# clean up table, look at AUC4 only
data <- data %>% select(TEVp, TEVs, rep, t1, t2, t3, t4, t5, AUC, min, sum)

# format data
data <- data %>%
  mutate(TEVs_pos = "P1'", # add column stating that P1` was mutated
         TEVs_AA = gsub("ENLYFQ", "", data$TEVs)) # write AA of mutated position in separate column)

# get data for some more TEVp variants (internal controls)
data_internal_controls <- read.table("./data/02_data_AUC_all.txt", sep = "\t", header = TRUE) %>% 
  # remove duplicate entries that stem from a mistake in previous data processing by Simon
  filter( !(TEVs == "ENLYFQS" & library.ID == "libLH016"),
          !(TEVs == "ENLYFQS" & library.ID == "libLH019")) %>%
  # only get the variants from the identical backbone plasmid
  filter(SphI.downstream.of.TEVs == "no") %>%
  # clean up
  select(rep, TEVp, TEVs, Note, min, t1, t2, t3, t4, t5, AUC4) %>% 
  rename(AUC = AUC4) 

# format internal control data likewise
data_internal_controls <- data_internal_controls %>%
  mutate(TEVs_pos = "P1'", # add column stating that P1` was mutated
         TEVs_AA = gsub("ENLYFQ", "", data_internal_controls$TEVs)) # write AA of mutated position in separate column)

# save data table
write.table(data, file = "./04_ABC_screen/results/data_clean.txt", sep = "\t", row.names=FALSE, quote = F)


# Summarize the data ------------------------------------------------------

  
# count how many unique TEVp variants and TEVp-TEVs combinations we have 
# including stop
unique_with_stop <- data %>%
  filter(min >= 20) %>% # only consider variants with >= threshold reads per time point
  group_by(rep) %>%
  summarize(n_TEVp_with_stop = n_distinct(TEVp),
            n_TEVp_TEVs_with_stop = n_distinct(paste(TEVp, TEVs, sep = "_")))

# count how many unique TEVp variants and TEVp-TEVs combinations we have 
# without stop

# Filter out TEVps containing stop codon
no_stop <- data %>%
  filter(min >= 20) %>% # only consider variants with >= threshold reads per time point
  filter(!grepl("\\*", TEVp))

unique_without_stop <- no_stop %>%
  group_by(rep) %>%
  summarize(n_TEVp_without_stop = n_distinct(TEVp),
            n_TEVp_TEVs_without_stop = n_distinct(paste(TEVp, TEVs, sep = "_")))


# Merge the two results based on 'rep'
summary <- merge(unique_with_stop, unique_without_stop, by = "rep", all.x = TRUE)


# save data table
write.table(summary, file = paste0("./04_ABC_screen/results/summary_min", 20, ".txt"), sep = "\t", row.names=FALSE, quote = F)


# make a data frame summarizing the coverage
coverage <- unique_without_stop %>%
  rename("TEVp library" = "n_TEVp_without_stop",
         "Combinatorial library" = "n_TEVp_TEVs_without_stop")

coverage <- pivot_longer(coverage, 
                         cols = c("TEVp library", "Combinatorial library"),
                         names_to = "library", 
                         values_to = "covered")

coverage$theoretical <- 0
coverage$theoretical[coverage$library == "TEVp library"] <- 8000
coverage$theoretical[coverage$library == "Combinatorial library"] <- 8000*20

# calculate some numbers
coverage <- coverage %>%
  mutate(missed = theoretical - covered,
         rel_coverage = covered / theoretical)

# pivot longer for plotting
coverage_longer <- coverage %>%
  select(rep, library, covered, missed) %>%
  pivot_longer(., cols = c(covered, missed), names_to = "category")

for (flask in c(1, 2)) {
  for (i in unique(coverage_longer$library)) {
    p <- ggplot(coverage_longer %>% filter(rep == flask,
                                           library == i), 
                aes(x = "", y = value, fill = category)) +
      geom_bar(stat = "identity", width = 1) +
      geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
      coord_polar(theta = "y") +
      labs(title = paste0(i, "\nReplicate ", flask, " min", 20)) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.title = element_blank()) +
      scale_fill_manual(values = c("missed" = "#cfcfcf", "covered" = "#69498a"))
    print(p)
    # save plot
    save_plot(p, filename = paste0("./04_ABC_screen/plots/pie_chart_min", 20, "_", i, "_", flask, "_.pdf"),
              base_height = 2, 
              base_asp = 1.3)
  }
}

## counting AAs ----

# do this for rep 1 only
flask <- 1

# get all the AAs
AAcount <- data %>%
  filter(rep == flask, 
         min >= 20) %>%
  select(TEVp, TEVs_AA) %>%
  mutate(Pos_30 = substr(TEVp, 1, 1),
         Pos_31 = substr(TEVp, 2, 2),
         Pos_32 = substr(TEVp, 3, 3)) %>%
  select(-TEVp)

# Reshape the data and count occurrences
count_s_df <- AAcount %>%
  pivot_longer(cols = names(AAcount), names_to = "Position", values_to = "AA") %>%
  group_by(AA, Position) %>%
  summarise(count = n())

# Load the RColorBrewer package
library(RColorBrewer)

# Create a custom color palette with 21 different colors
my_palette <- brewer.pal(10, "Paired")
my_colors <- c(my_palette, my_palette, my_palette[1])

# Create the stacked bar plot with the custom color palette
p1 <- ggplot(count_s_df, aes(x = Position, y = count, fill = AA)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +  # Set the custom color palette
  labs(x = "Position", y = "Count", fill = NULL) +
  theme_cowplot(font_size = 10) +
  theme(legend.position = "right") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_x_discrete(labels = c("30", "31", "32", "P1'")) +
  theme(legend.text = element_text(size = 10))

# Modify the y-axis labels to display 'k' for 1000
p1 <- p1 + scale_y_continuous(labels = scales::label_number(scale = 1e-3, accuracy = 1, suffix = "k"))

print(p1)

## plot theoretical AA distribution ----

# Load necessary library
library(Biostrings)

# Define the bases for N and K
N <- c("A", "T", "C", "G")
K <- c("G", "T")

# Generate all possible NNK codons
NNK_codons <- c()
for (base1 in N) {
  for (base2 in N) {
    for (base3 in K) {
      NNK_codons <- c(NNK_codons, paste0(base1, base2, base3))
    }
  }
}

# Translate the NNK codons to amino acids
amino_acids <- translate(DNAStringSet(NNK_codons))

# Calculate the frequency of each amino acid
amino_acid_freq <- table(amino_acids) / length(NNK_codons)

# Convert amino_acid_freq to a data frame
amino_acid_freq_df <- as.data.frame(amino_acid_freq) 
colnames(amino_acid_freq_df) <- c("AA", "Frequency")


p2 <- ggplot(amino_acid_freq_df, aes(x = "NNK", y = Frequency, fill = AA)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +  # Set the custom color palette
  labs(y = "Theoretical frequency", fill = NULL) +
  theme_cowplot(font_size = 10) +
  theme(legend.position = "right") +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.text = element_text(size = 10))


# Put all plots into a grid
pg <- plot_grid(p2 + theme(legend.position = "none"), p1, rel_widths = c(1,3))
pg

# save grid
save_plot(pg, filename = "./04_ABC_screen/plots/AA_count_vs_theory.pdf", 
          base_height = 3.71,
          base_asp = 1.4) 


# 2_noise correction ------------------------------------------------------

# get information about noise
df_noise <- read.table("./01_internal_controls/results/df_noise.txt", sep = "\t", header = TRUE)

cutoff <- filter(df_noise, outlier_removed == F) %>% pull(cutoff)
mu <- filter(df_noise, outlier_removed == F) %>% pull(mu)

## manipulate data acc. to noise ##

# calculate noise-blanked values (delta noise AUC, dnAUC)
# for that, set all values below the cutoff to the mean of the noise
# then blank to the mean of the noise
data$dnAUC <- ifelse(data$AUC <= cutoff, mu, data$AUC)-mu

# calculate noise-blanked values rel to wt TEVs (dnAUCrel)
# get values of noise-blanked wt TEVs
nblanked_ENLYFQS <- data %>% filter(TEVs == "ENLYFQS") %>% select(rep, TEVp, dnAUC) %>% rename("nblanked_ENLYFQS" = "dnAUC")
# join data with these values
data <- full_join(data, nblanked_ENLYFQS, by = c("TEVp", "rep")) 
# calculate normalized blanked values
data$dnAUCrel <- data$dnAUC / data$nblanked_ENLYFQS


# 3_identification of interesting variants ----------------------------------

## orthogonal variants ----

# to identify orthogonal variants divide the first maximum by the second maximum
# within each TEVp
orthogonality <- data %>%
  filter(min >= 20) %>%
  group_by(TEVp, rep) %>%
  arrange(desc(AUC)) %>%
  summarize(max1_AUC = max(AUC),
            max2_AUC = nth(AUC, 2),
            orthog = max1_AUC / max2_AUC,
            top_TEVs_AA = TEVs_AA[which.max(AUC)],
            substrates = n()) %>%
  arrange(desc(orthog)) %>%
  # remove variants containing a stop codon
  filter(!(grepl("\\*", TEVp)))

# save data table
write.table(orthogonality, file = "./04_ABC_screen/results/orthogonality.txt", sep = "\t", row.names=FALSE, quote = F)

## promiscuous variants ----

# install and load the entropy package 
library(entropy)

# to identify promiscuous variants calculate the entropy of AUC within each TEVp
# and multiply it by the mean activity of that TEVp
promiscuity <- data %>%
  filter(min >= 20) %>%
  group_by(TEVp, rep) %>%
  summarize(promisc = entropy(AUC) * mean(AUC),
            substrates = n()) %>%
  arrange(desc(promisc)) %>%
  # remove variants containing a stop codon
  filter(!(grepl("\\*", TEVp)))


# save data table
write.table(promiscuity, file = "./04_ABC_screen/results/promiscuity.txt", sep = "\t", row.names=FALSE, quote = F)

# 4_heat maps ---------------------------------------------------------------

## to plot all variants ----

# write function to generate P1' heat maps
heatmap_all <- function(
  replicate, TEVp_var, onlyactive = FALSE, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data,
                      min >= 20,
                      rep == replicate,
                      TEVp %in% TEVp_var) %>%
    # manually remove the one extreme outlier
    filter(!(TEVp == "VDS" & TEVs == "ENLYFQT"))
  
  # add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing values
  # create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations - 3381 lines
  df_combinations <- expand.grid(
    TEVp = unique(plot_data$TEVp),
    TEVs_AA = unique(plot_data$TEVs_AA))
  # merge with data of actually measured TEVp x TEVs combinations
  plot_data <- plot_data %>%
    full_join(., df_combinations,
              by = c("TEVp", "TEVs_AA")) %>%
    # give the name "values" to the column that should be plotted and is defined by "activity"
    rename(values = all_of(activity))
  
  # print max value
  cat("Maximum value in this plot:", max(plot_data$values, na.rm = TRUE))
  
  # manually define the order of TEVp variants
  TEVp_levels <- plot_data %>%
    group_by(TEVp) %>% 
    summarize(activity_mean = mean(na.omit(values))) %>% 
    # slice_sample(., prop = 1) %>% # Shuffle rows randomly
    {if (onlyactive) filter(., activity_mean > 0) else .} %>% # conditional filter for active variants
    arrange(desc(activity_mean)) %>% 
    pull(TEVp)
  
  # plotting
  p <- ggplot(plot_data, aes(x = factor(TEVp, levels = TEVp_levels), 
                             y = factor(TEVs_AA, levels = rev(sort(unique(TEVs_AA)))),
                             fill = values)) +
    geom_tile() + # make heat map
    scale_fill_viridis(option="magma",
                       # limits = c(0, 0.225) # set a manual limit to be consistent with other plots in the figure
    ) +
    theme_cowplot(font_size = 12) +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +  # Remove x-axis tick marks
    labs(fill = legend) +
    scale_y_discrete(name = "Amino acid in P1'", expand = c(0, 0)) +
    # scale_y_reverse() + # reverse the y-axis
    scale_x_discrete(name = paste0("TEVp variants (n = ",
                                   format(length(TEVp_levels), big.mark = ","),      
                                   ")"), expand = c(0, 0),
                     # breaks = TEVp_levels[seq(0, length(TEVp_levels), 500)]
                     )
  
  # show plot
  print(p)
  
  # save plot to file
  save_plot(p, filename = paste0("./04_ABC_screen/plots/heatmaps/", file_name, ".png"),
            base_height = height,
            base_asp = asp_ratio)
}



## plotting of all variants ----

# heat map of all variants showing noise-blanked AUC values
heatmap_all(replicate = 1,
            TEVp_var = data$TEVp[!grepl("\\*", data$TEVp)] %>% unique(), # remove variants containing a stop codon
            activity = "dnAUC",
            legend = "Activity",
            file_name = "screen_all_variants_dnAUC", 
            height = 3.15, 
            asp_ratio = 1.2)

# heat map of variants active on at least one substrate showing noise-blanked AUC values
heatmap_all(replicate = 1,
            TEVp_var = data$TEVp[!grepl("\\*", data$TEVp)] %>% unique(), # remove variants containing a stop codon
            onlyactive = TRUE,
            activity = "dnAUC",
            legend = "Activity",
            file_name = "screen_active_variants_dnAUC", 
            height = 3.15, 
            asp_ratio = 1.2)



# heat map of variants with at least 15 TEVss showing noise-blanked AUC values
heatmap_all(replicate = 1,
            TEVp_var = data %>% 
              filter(rep == 1) %>% 
              group_by(TEVp) %>% 
              summarize(TEVpxTEVs = n()) %>% 
              filter(TEVpxTEVs >=15) %>% 
              pull(TEVp),
            activity = "dnAUC",
            legend = "dnAUC",
            file_name = "screen_15substrates_dnAUC", 
            height = 3.15, 
            asp_ratio = 3)

## to plot potential hits ----

# write function to generate P1' heat maps of promiscuous hits
heatmap_promiscuous <- function(
  replicate, TEVp_var, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data,
                      min >= 20,
                      rep == replicate,
                      TEVp %in% TEVp_var)
  
  # add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing values
  # create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations - 3381 lines
  df_combinations <- expand.grid(
    TEVp = unique(plot_data$TEVp),
    TEVs_AA = unique(plot_data$TEVs_AA))
  # merge with data of actually measured TEVp x TEVs combinations
  plot_data <- plot_data %>%
    full_join(., df_combinations,
              by = c("TEVp", "TEVs_AA")) %>%
    # give the name "values" to the column that should be plotted and is defined by "activity"
    rename(values = all_of(activity))
  
  # manually define the order of TEVp variants
  TEVp_levels <- plot_data %>%
    group_by(TEVp) %>% 
    summarize(activity_mean = mean(na.omit(values))) %>% 
    arrange(desc(activity_mean)) %>% 
    pull(TEVp)
  
  # plotting
  p <- ggplot(plot_data, aes(x = factor(TEVp, levels = TEVp_levels), 
                             y = factor(TEVs_AA, levels = rev(sort(unique(TEVs_AA)))),
                             fill = values)) +
    geom_tile() + # make heat map
    scale_fill_viridis(option="magma") +
    theme_cowplot(font_size = 12) +
    theme(axis.line = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("") +
    labs(fill = legend) +
    scale_y_discrete(name = "Amino acid in P1'", expand = c(0, 0))
  
  # show plot
  print(p)
  
  # save plot to file
  save_plot(p, filename = paste0("./04_ABC_screen/plots/heatmaps/promiscuous/", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}


# write function to generate P1' heat maps of orthogonal hits
heatmap_orthogonal <- function(
  replicate, TEVp_var, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data,
                      min >= 20,
                      rep == replicate,
                      TEVp %in% TEVp_var)
  
  # add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing values
  # create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations - 3381 lines
  df_combinations <- expand.grid(
    TEVp = unique(plot_data$TEVp),
    TEVs_AA = unique(plot_data$TEVs_AA))
  # merge with data of actually measured TEVp x TEVs combinations
  plot_data <- plot_data %>%
    full_join(., df_combinations,
              by = c("TEVp", "TEVs_AA")) %>%
    # give the name "values" to the column that should be plotted and is defined by "activity"
    rename(values = all_of(activity))
  
  # manually define the order of TEVp variants, order acc. to activity of top TEVs_AA
  TEVp_levels <- plot_data %>%
    arrange(desc(values)) %>% 
    pull(TEVp) %>% unique()
  
  # plotting
  p <- ggplot(plot_data, aes(x = factor(TEVp, levels = TEVp_levels), 
                             y = factor(TEVs_AA, levels = rev(sort(unique(TEVs_AA)))),
                             fill = values)) +
    geom_tile() + # make heat map
    scale_fill_viridis(option="magma", limits = c(0, 0.25)) +
    theme_cowplot(font_size = 12) +
    theme(axis.line = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("") +
    labs(fill = legend) +
    scale_y_discrete(name = "Amino acid in P1'", expand = c(0, 0))
  
  # show plot
  print(p)
  
  # save plot to file
  save_plot(p, filename = paste0("./04_ABC_screen/plots/heatmaps/orthogonal/", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}


## plotting of potential hits ----

# heat map of top promiscuous variants for each library
heatmap_promiscuous(replicate = 1,
                    TEVp_var = promiscuity %>%
                      filter(rep == 1,
                             substrates >= 15) %>%
                      arrange(desc(promisc)) %>%
                      head(50) %>% 
                      pull(TEVp),
                    activity = "dnAUC",
                    legend = "dnAUC",
                    file_name = "15substrates_dnAUC", 
                    height = 4.15, 
                    asp_ratio = 3)


# heat map of top orthogonal variants
heatmap_orthogonal(replicate = 1,
                   TEVp_var = orthogonality %>%
                     filter(rep == 1,
                            substrates >= 15) %>%
                     arrange(desc(orthog)) %>%
                     head(50) %>% 
                     pull(TEVp),
                   activity = "dnAUC",
                   legend = "dnAUC",
                   file_name = "15substrates_dnAUC_best", 
                   height = 4.15, 
                   asp_ratio = 2.5)


# heat map of top orthogonal variants for each TEVs
# define the TEVs_AAs
TEVs_AA <- data %>% pull(TEVs_AA) %>% unique()

# create folders for the libraries
dir.create("./04_ABC_screen/plots/heatmaps/orthogonal/per_substrate/")

# create heat maps  
for (AA in TEVs_AA) {
  replic <- 1
  # get TEVps
  TEVps <- orthogonality %>%
    filter(rep == replic,
           substrates >= 15,
           top_TEVs_AA == AA) %>% #,
           # max1_AUC >= 0.1) %>% # set a threshold for the activity
    arrange(desc(orthog)) %>%
    head(50) %>% 
    pull(TEVp)
  # Check if there are no entries for the current TEVp, in that case skip this current variable
  if (length(TEVps) == 0) {
    next  # Skip the current iteration and move to the next variable (this prevents the loop from stopping)
  }
  # use function to make plot
  heatmap_orthogonal(replicate = replic,
                     TEVp_var = TEVps,
                     activity = "dnAUC",
                     legend = "dnAUC",
                     file_name = paste0("per_substrate/15substrates_dnAUC_", AA), 
                     height = 4.15, 
                     asp_ratio = 2.5)
}


## to plot selected hits ----

# write function to generate P1' heat maps of manually selected potential hits
# i.e. of variants that seem especially orthogonal or promiscuous
heatmap_hits <- function(
  replicate, TEVp_var, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data,
                      min >= 10,
                      rep == replicate,
                      TEVp %in% TEVp_var)
  # add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing values
  # create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations - 3381 lines
  df_combinations <- expand.grid(
    TEVp = unique(plot_data$TEVp),
    TEVs_AA = unique(plot_data$TEVs_AA))
  # merge with data of actually measured TEVp x TEVs combinations
  plot_data <- plot_data %>%
    full_join(., df_combinations,
              by = c("TEVp", "TEVs_AA")) %>%
    # give the name "values" to the column that should be plotted and is defined by "activity"
    rename(values = all_of(activity))
  
  # manually specify a category for each TEVp
  plot_data$category <- "specific"
  plot_data$category[plot_data$TEVp == 'TSL'] <- "wt"
  plot_data$category[plot_data$TEVp %in%  c('GST', 'QSQ', 'VST')] <- "subset"
  plot_data$category[plot_data$TEVp %in%  c('ASL', 'ATL', 'DSL')] <- "promisc."
  
  # manually define the order of TEVp variants
  TEVp_levels <- c('DSL', 'ATL', 'ASL', 'TSL', 'VST', 'QSQ', 'GST', 'GGT', 'PQW', 
                   'RPD', 'DRA', 'QRS', 'ERA', 'IGR', 'HGR', 'RAV', 'HPG', 'AGM',
                   'LAH', 'RWG', 'CGA', 'TCR')
  
  # change facet order
  plot_data$category <- factor(plot_data$category, levels = c("specific", "subset", "wt", "promisc."))
  
  # plotting
  p <- ggplot(plot_data, aes(x = TEVs_AA, 
                             y = factor(TEVp, levels = TEVp_levels), 
                             fill = values)) +
    geom_tile() + # make heat map
    scale_fill_viridis(option="magma") +
    theme_cowplot(font_size = 12) +
    # theme(legend.position = "top") +
    theme(axis.line = element_blank()) +
    theme(axis.title.y = element_blank()) +
    #    theme(axis.text.y = element_text(family = "Courier")) + # use a monospace font for the y-axis labels
    labs(fill = legend) +
    scale_x_discrete(name = "Amino acid in P1'", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    #    theme(axis.text.x = element_text(size = 10)) + # change size of x axis labels
    facet_grid(category ~ ., scales = "free", space = "free") # switch rows and columns
  
  # show plot
  print(p)
  # save plot to file
  save_plot(p, filename = paste0("./04_ABC_screen/plots/heatmaps/", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}


## plotting of selected hits ----

# specify the manually identified hits
hits <- c('DSL', 'ATL', 'ASL', 'TSL', 'VST', 'QSQ', 'GST', 'GGT', 'PQW', 
          'RPD', 'DRA', 'QRS', 'ERA', 'IGR', 'HGR', 'RAV', 'HPG', 'AGM', 
          'LAH', 'RWG', 'CGA', 'TCR')

# heat map of all variants showing noise corrected AUC values
heatmap_hits(replicate = "1",
             TEVp_var = hits,
             activity = "dnAUC",
             legend = "dnAUC",
             file_name = paste0("hits_dnAUC_rep1"),
             height = 5, 
             asp_ratio = 1)







