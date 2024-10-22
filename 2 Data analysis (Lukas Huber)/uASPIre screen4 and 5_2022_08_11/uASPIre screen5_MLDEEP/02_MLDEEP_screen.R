# load libraries
library(tidyverse)
library(viridis)
library(cowplot)

# create folders ----------------------------------------------------------

# create a folder to store the individual plots and data in
dir.create("./02_MLDEEP")
dir.create("./02_MLDEEP/results")
dir.create("./02_MLDEEP/plots")
dir.create("./02_MLDEEP/plots/heatmaps")
dir.create("./02_MLDEEP/plots/heatmaps/orthogonal")
dir.create("./02_MLDEEP/plots/heatmaps/promiscuous")

# prepare data ------------------------------------------------------------

# read data
data <- read.table("./data/Data_MLD.txt", sep = "\t", header = TRUE)

# clean up table
data <- data %>% 
  select(-source_TEVp, -source_TEVs)

# format data
data <- data %>%
  mutate(TEVs_pos = "P1'", # add column stating that P1` was mutated
         TEVs_AA = gsub("ENLYFQ", "", data$TEVs), # write AA of mutated position in separate column)
         TEVp = paste0(TEVp_ABC, "_", TEVp_DEF, "_", TEVp_XYZ)) 

# remove TEVp variants with duplicate entries that stem from the way how Simon assigned the library to variants
# In Simon's code, a variant is assigned to a library according to its DNA code, not AA code
# Thus, synonymous mutants yield duplicate entries

# check for duplicate TEVps
check <- data %>%
  group_by(TEVp, TEVs_AA) %>%
  summarize(n = n()) %>%
  filter(n >= 2) %>%
  arrange(desc(n))

# get the data of duplicated TEVps
duplicates <- data %>% filter(TEVp %in% check$TEVp)

# filter for entries that are assigned to one of the six-position libraries 
to_be_removed <- duplicates %>% filter(library %in% c("ABCXYZ", "ABCDEF"))

# remove those entries from the data
data <- anti_join(data, to_be_removed)

# get data for some more TEVp variants (internal controls)
data_internal_controls <- read.table("./data/Data_MLD_internal_controls.txt", sep = "\t", header = TRUE) %>%  # clean up
  select(TEVp, TEVs, Note, min, t1, t2, t3, t4, t5, AUC)

# format internal control data likewise
data_internal_controls <- data_internal_controls %>%
  mutate(TEVs_pos = "P1'", # add column stating that P1` was mutated
         TEVs_AA = gsub("ENLYFQ", "", data_internal_controls$TEVs)) # write AA of mutated position in separate column)


## dealing with noise (dnAUC, dnAUCrel)----

# get information about noise
df_noise <- read.table("./01_internal_controls/results/df_noise.txt", sep = "\t", header = TRUE)

cutoff <- df_noise %>% pull(cutoff)
mu <- df_noise %>% pull(mu)

## manipulate data acc. to noise ##

# calculate noise-blanked values (delta noise AUC, dnAUC)
# for that, set all values below the cutoff to the mean of the noise
# then blank to the mean of the noise
data$dnAUC <- ifelse(data$AUC <= cutoff, mu, data$AUC)-mu

# calculate noise-blanked values rel to wt TEVs (dnAUCrel)
# get values of noise-blanked wt TEVs
nblanked_ENLYFQS <- data %>% filter(TEVs == "ENLYFQS") %>% select(TEVp, dnAUC) %>% rename("nblanked_ENLYFQS" = "dnAUC")
# join data with these values
data <- full_join(data, nblanked_ENLYFQS, by = c("TEVp")) 
# calculate normalized blanked values
data$dnAUCrel <- data$dnAUC / data$nblanked_ENLYFQS

# remove samples above the positive control as they might very likely be false positives
posctrl <- data_internal_controls %>%
  filter(TEVp == "0 C151A no ssrA") %>% 
  pull(AUC)

data <- data %>% filter(AUC <= posctrl)

# save data table
write.table(data, file = "./02_MLDEEP/results/data_clean.txt", sep = "\t", row.names=FALSE, quote = F)

# Summarize the data ------------------------------------------------------

# define a read threshold
threshold <- 20

# count how many unique TEVp variants and TEVp-TEVs combinations we have 
# including stop
unique_with_stop <- data %>%
  rename(space_library = library) %>%
  filter(space_library != "WT") %>%
  filter(min >= threshold) %>% # only consider variants with >= threshold reads per time point
  group_by(space_library) %>%
  summarize(n_TEVp_with_stop = n_distinct(TEVp),
            n_TEVp_TEVs_with_stop = n_distinct(paste(TEVp, TEVs, sep = "_")))

# count how many unique TEVp variants and TEVp-TEVs combinations we have 
# without stop

# Filter out TEVps containing stop codon
no_stop <- data %>%
  rename(space_library = library) %>%
  filter(space_library != "WT") %>%
  filter(min >= threshold) %>% # only consider variants with >= threshold reads per time point
  filter(!grepl("\\*", TEVp))

unique_without_stop <- no_stop %>%
  filter(space_library != "WT") %>%
  group_by(space_library) %>%
  summarize(n_TEVp_without_stop = n_distinct(TEVp),
            n_TEVp_TEVs_without_stop = n_distinct(paste(TEVp, TEVs, sep = "_")))


# Merge the two results based on 'space_library'
summary <- merge(unique_with_stop, unique_without_stop, by = "space_library", all.x = TRUE)

# save data table
write.table(summary, file = paste0("./02_MLDEEP/results/summary_min", threshold, ".txt"), sep = "\t", row.names=FALSE, quote = F)


# make a data frame summarizing the coverage
coverage <- unique_without_stop %>%
  rename("TEVp library" = "n_TEVp_without_stop",
         "Combinatorial library" = "n_TEVp_TEVs_without_stop")

coverage <- pivot_longer(coverage, 
                         cols = c("TEVp library", "Combinatorial library"),
                         names_to = "library", 
                         values_to = "covered")

# add theoretical library sizes
coverage$theoretical <- 0
coverage$theoretical[coverage$library == "TEVp library" & 
                       nchar(coverage$space_library) == 3] <- 20^3
coverage$theoretical[coverage$library == "Combinatorial library" & 
                       nchar(coverage$space_library) == 3] <- 20^3*20
coverage$theoretical[coverage$library == "TEVp library" & 
                       nchar(coverage$space_library) == 6] <- 20^6
coverage$theoretical[coverage$library == "Combinatorial library" & 
                       nchar(coverage$space_library) == 6] <- 20^6*20
coverage$theoretical[coverage$library == "TEVp library" & 
                       nchar(coverage$space_library) == 9] <- 20^9
coverage$theoretical[coverage$library == "Combinatorial library" & 
                       nchar(coverage$space_library) == 9] <- 20^9*20

# calculate some numbers
coverage <- coverage %>%
  mutate(missed = theoretical - covered,
         rel_coverage = covered / theoretical)

# pivot longer for plotting
coverage_longer <- coverage %>%
  select(space_library, library, covered, missed) %>%
  pivot_longer(., cols = c(covered, missed), names_to = "category")

for (space in unique(coverage_longer$space_library)) {
  for (i in unique(coverage_longer$library)) {
    p <- ggplot(coverage_longer %>% filter(space_library == space,
                                           library == i), 
                aes(x = "", y = value, fill = category)) +
      geom_bar(stat = "identity", width = 1) +
      geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
      coord_polar(theta = "y") +
      labs(title = paste0(i, "\nSpace ", space, " min", threshold)) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.title = element_blank()) +
      scale_fill_manual(values = c("missed" = "#cfcfcf", "covered" = "#69498a"))
    print(p)
    # save plot
    save_plot(p, filename = paste0("./02_MLDEEP/plots/pie_chart_min", threshold, "_", i, "_", space, ".pdf"),
              base_height = 2, 
              base_asp = 1.3)
  }
}




# 1_identification of interesting variants ----------------------------------

## orthogonal variants ----

# to identify orthogonal variants divide the first maximum by the second maximum
# within each TEVp
orthogonality <- data %>%
  group_by(TEVp, library) %>%
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
write.table(orthogonality, file = "./02_MLDEEP/results/orthogonality.txt", sep = "\t", row.names=FALSE, quote = F)
  
## promiscuous variants ----

# install and load the entropy package 
library(entropy)

# to identify promiscuous variants calculate the entropy of AUC within each TEVp
# and multiply it by the mean activity of that TEVp
promiscuity <- data %>%
  group_by(TEVp, library) %>%
  summarize(promisc = entropy(AUC) * mean(AUC),
            substrates = n()) %>%
  arrange(desc(promisc)) %>%
  # remove variants containing a stop codon
  filter(!(grepl("\\*", TEVp)))

# save data table
write.table(promiscuity, file = "./02_MLDEEP/results/promiscuity.txt", sep = "\t", row.names=FALSE, quote = F)

# 2_heat maps ---------------------------------------------------------------

## to plot heat maps generally ----

# write function to generate P1' heat maps
heatmap_all <- function(
  TEVp_var, onlyactive = FALSE, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data,
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
  
  # print max value
  cat("Maximum value in this plot:", max(plot_data$values, na.rm = TRUE))
  
  # manually define the order of TEVp variants
  TEVp_levels <- plot_data %>%
    group_by(TEVp) %>% 
    summarize(activity_mean = mean(na.omit(values))) %>% 
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
  save_plot(p, filename = paste0("./02_MLDEEP/plots/heatmaps/", file_name, ".png"),
            base_height = height,
            base_asp = asp_ratio)
}



## plotting of all variants ----

# heat map of all variants showing noise-blanked AUC values
# for each library individually

# # pull out all different libraries
# libraries <- data %>% pull(library) %>% unique() %>% .[! . == "WT"]

# define libraries for which we want to generate plots
libraries <- c("XYZ", "ABCXYZ")

# do this for each library
for (lib in libraries) {
  heatmap_all(TEVp_var = data %>%
                filter(library == lib,
                       !(grepl("\\*", TEVp))) %>% # remove variants containing a stop codon
                pull(TEVp) %>%
                unique(), 
              activity = "dnAUC",
              legend = "Activity",
              file_name = paste0("all_variants_dnAUC_", lib), 
              height = 3.15, 
              asp_ratio = 1.2)
}

# do this for each library, plotting only variants active on at least one substrate
for (lib in libraries) {
  heatmap_all(TEVp_var = data %>%
                filter(library == lib,
                       !(grepl("\\*", TEVp))) %>% # remove variants containing a stop codon
                pull(TEVp) %>%
                unique(), 
              onlyactive = TRUE,
              activity = "dnAUC",
              legend = "Activity",
              file_name = paste0("active_variants_dnAUC_", lib), 
              height = 3.15, 
              asp_ratio = 1.2)
}


## to plot heat maps of hits ----

# write function to generate P1' heat maps of promiscuous hits
heatmap_promiscuous <- function(
  TEVp_var, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data,
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
  save_plot(p, filename = paste0("./02_MLDEEP/plots/heatmaps/", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}


# write function to generate P1' heat maps of orthogonal hits
heatmap_orthogonal <- function(
  TEVp_var, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data,
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
  save_plot(p, filename = paste0("./02_MLDEEP/plots/heatmaps/", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}

## plotting of potential hits ----

# heat map of top promiscuous variants for each library
for (lib in libraries) {
  heatmap_promiscuous(TEVp_var = promiscuity %>%
                        filter(library == lib,
                               substrates >= 15) %>%
                        arrange(desc(promisc)) %>%
                        head(50) %>% 
                        pull(TEVp),
                      activity = "dnAUC",
                      legend = "dnAUC",
                      file_name = paste0("promiscuous/15substrates_dnAUC_", lib), 
                      height = 4.15, 
                      asp_ratio = 2.5)
}


# heat map of top orthogonal variants with high activity for each library 
for (lib in libraries) {
  heatmap_orthogonal(TEVp_var = orthogonality %>%
                       filter(library == lib,
                              substrates >= 15,
                              max1_AUC >= 0.1) %>%
                       arrange(desc(orthog)) %>%
                       head(50) %>% 
                       pull(TEVp),
                     activity = "dnAUC",
                     legend = "dnAUC",
                     file_name = paste0("orthogonal/15substrates_dnAUC_active_", lib), 
                     height = 4.15, 
                     asp_ratio = 2.5)
}

