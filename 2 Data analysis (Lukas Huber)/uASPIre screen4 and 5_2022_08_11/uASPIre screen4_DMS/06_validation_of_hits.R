### THIS SCRIPT IS ABOUT THE VALIDATION OR RETESTING OF VARIANTS THAT HAVE BEEN
### IDENTIFIED IN THE ABC SCREEN AS PUTATIVELY PROMISCUOUS OR SPECIFIC FOR A 
### CERTAIN AMINO ACID IN P1'. THOSE VARIANTS HAVE BEEN CLONED AGAIN AND TESTED
### AGAIN. SOME OF THE TEVp-TEVs COMBINATIONS HAVE ALSO BEEN TESTED IN VITRO. 
### AT THE END, THE IN VITRO RESULTS ARE ALSO COMPARED TO THE RECORDER RESULTS. 

# load libraries
library(tidyverse)
library(viridis)
library(cowplot)
library(ggseqlogo)
library(RColorBrewer)



# create folders ----------------------------------------------------------

# create a folder to store the individual plots in
dir.create("./06_validation_of_hits")
dir.create("./06_validation_of_hits/plots")
dir.create("./06_validation_of_hits/plots/heatmaps")

# prepare data ------------------------------------------------------------

## full specificity (using PacBio look up table) ----

# read data
data <- read.table("./data/Data_DMS_validation_PacBio.txt", sep = "\t", header = TRUE)

## do the following to format the table
# select only single mutants of TEVs
TEVs_singlemut <- data %>% filter(n_mut_TEVs == 1) %>%
  select(TEVp, TEVs, t1, t2, t3, t4, t5, AUC)

# define which position in TEVs is mutated
TEVs_singlemut <- TEVs_singlemut %>% mutate(TEVs_pos = mapply(function(x, y) which(x != y)[1], 
       strsplit(TEVs_singlemut$TEVs, ""), strsplit("ENLYFQS", "")))

# write AA of mutated position in separate column
TEVs_singlemut <- TEVs_singlemut %>% mutate(TEVs_AA = substring(TEVs, TEVs_pos, TEVs_pos))

# rename the positions using the Schechter and Berger nomenclature
TEVs_singlemut$TEVs_pos[TEVs_singlemut$TEVs_pos == '1'] <- "P6"
TEVs_singlemut$TEVs_pos[TEVs_singlemut$TEVs_pos == '2'] <- "P5"
TEVs_singlemut$TEVs_pos[TEVs_singlemut$TEVs_pos == '3'] <- "P4"
TEVs_singlemut$TEVs_pos[TEVs_singlemut$TEVs_pos == '4'] <- "P3"
TEVs_singlemut$TEVs_pos[TEVs_singlemut$TEVs_pos == '5'] <- "P2"
TEVs_singlemut$TEVs_pos[TEVs_singlemut$TEVs_pos == '6'] <- "P1"
TEVs_singlemut$TEVs_pos[TEVs_singlemut$TEVs_pos == '7'] <- "P1'"

# add the two new columns to the original data set
data <- full_join(data, TEVs_singlemut)

# clean up data set
data <- data %>% select(TEVp, TEVs, n_mut_TEVs, TEVs_pos, TEVs_AA, t1, t2, t3, t4, t5, sum, min, AUC)


## P1' specificity (using manual look up table, i.e. of 20 manually cloned P1' backbones) ----

# read data
data_manual <- read.table("./data/Data_DMS_validation_manual.txt", sep = "\t", header = TRUE)

# format data 
data_manual <- data_manual %>% 
  mutate(TEVp_parent = substring(TEVp, 1, 1), # separate the name of TEVp into parent and position 30-32
         TEVp = substring(TEVp, 3, 5))

# get data for some more TEVp variants (internal controls)
data_internal_controls <- read.table("./data/Data_DMS_internal_controls.txt", sep = "\t", header = TRUE)

# add parent variants to the data
data_manual <- data_internal_controls %>%
  select(-Note) %>%
  filter(., TEVp %in% c("0", "I")) %>% # get parent variants TEVp 0 and I
  rename("TEVp_parent" = "TEVp") %>% 
  mutate(TEVp = "TSL")  %>% # add information about AAs in pos 30-32
  rbind(., data_manual) # bind with other data

# format data
data_manual <- data_manual %>% 
  mutate(TEVs_pos = "P1'", # add column stating that P1` was mutated
         TEVs_AA = gsub("ENLYFQ", "", data_manual$TEVs)) # write AA of mutated position in separate column)

# clean up data set 
data_manual <- data_manual %>% select(TEVp_parent, TEVp, TEVs, TEVs_pos, TEVs_AA, t1, t2, t3, t4, t5, sum, min, AUC)


### correction for noise ----

# get information about noise of this screen
noise <- read.table("./01_internal_controls/results/df_noise.txt", sep = "\t", header = TRUE)
cutoff <- noise$cutoff
mu <- noise$mu

# calculate noise-blanked values (delta noise AUC)
# for that, set all values below the cutoff to the mean of the noise
# then blank to the mean of the noise
data_manual$dnAUC <- ifelse(data_manual$AUC <= cutoff, mu, data_manual$AUC)-mu

### noise blanked values rel to wt TEVs
# get values of noise-blanked wt TEVs
nblanked_ENLYFQS <- data_manual %>% filter(TEVs == "ENLYFQS") %>% select(TEVp_parent, TEVp, dnAUC) %>% rename("nblanked_ENLYFQS" = "dnAUC")
# join data with these values
data_manual <- full_join(data_manual, nblanked_ENLYFQS, by = c("TEVp_parent", "TEVp")) 
# calculate normalized blanked values
data_manual <- data_manual %>% mutate(dnAUCrel = dnAUC / nblanked_ENLYFQS)


# Identify duplicates based on all columns
duplicates <- duplicated(data_manual)

# View the duplicated rows
duplicated <- data_manual[duplicates,]

# prepare plotting --------------------------------------------------------


## to plot heat maps of full specificity profiles ----

# manually define the order of TEVp variants
TEVp_levels <- data$TEVp %>% unique()

# write function to generate heat maps
heatmap <- function(
  TEVp_var, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- rbind(
    # first, copy data of wildtype TEVs once for each TEVs position and append it
    # this step is necessary to avoid "holes" in heat maps
    filter(data, n_mut_TEVs != 0),
    filter(data, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P6", TEVs_AA = "E"),
    filter(data, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P5", TEVs_AA = "N"),
    filter(data, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P4", TEVs_AA = "L"),
    filter(data, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P3", TEVs_AA = "Y"),
    filter(data, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P2", TEVs_AA = "F"),
    filter(data, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P1", TEVs_AA = "Q"),
    filter(data, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P1'", TEVs_AA = "S")) %>%
    # filter for TEVp variants
    filter(n_mut_TEVs <= 1, TEVs_AA != "*", TEVp %in% TEVp_var)
  
  # add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing values
  # create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations - 3381 lines
    df_combinations <- expand.grid(
      TEVp = unique(plot_data$TEVp),
      TEVs_pos = unique(plot_data$TEVs_pos),
      TEVs_AA = unique(plot_data$TEVs_AA))
    # merge with data of actually measured TEVp x TEVs combinations
    plot_data <- plot_data %>%
      full_join(., df_combinations,
                by = c("TEVp", "TEVs_pos", "TEVs_AA")) %>%
      # make the positions as factors and define their order for the heat map
      mutate(TEVs_pos = factor(TEVs_pos, levels=c("P1'", "P1", "P2", "P3", "P4", "P5", "P6"))) %>%
      # give the name "values" to the column that should be plotted and is defined by "activity"
      rename(values = all_of(activity))
    
  # plotting
  p <- ggplot(plot_data, aes(x = TEVs_AA, y = TEVs_pos, fill = values)) +
    geom_tile() +
    facet_wrap(~factor(TEVp, levels = TEVp_levels)) +
    scale_fill_viridis(option="magma") +
    theme_cowplot(font_size = 10) +
    theme(axis.line = element_blank()) +
    theme(axis.title.y = element_blank()) +
    labs(fill = legend) +
    xlab("Amino acid") +
    theme(axis.text.x = element_text(size = 8)) # change size of x axis labels
  # show plot
  print(p)
  # save plot to file
  save_plot(p, filename = paste0("./06_validation_of_hits/plots/heatmaps/", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}


## to plot heat maps of P1' specificity profiles ----

# write function to generate P1' heat maps
heatmap_P1p <- function(
  TEVp_var, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data_manual,
                      TEVp_parent == "I",
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
  save_plot(p, filename = paste0("./06_validation_of_hits/plots/heatmaps/P1p_", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}


# Full specificity heat maps (PacBio) ------------------------------------------

## raw AUCs ----

# heat map of all variants showing raw AUC values
heatmap(TEVp_var = pull(data, TEVp),
        activity = "AUC",
        legend = "AUC",
        file_name = "all_variants_AUC",
        height = 3.71*2,
        asp_ratio = 1.618)


# P1' heat maps (manual) -------------------------------------------------------

hits <- data_manual %>% filter(TEVp != "ACR") %>% pull(TEVp) %>% unique()

# heat map of all variants showing raw AUC values
heatmap_P1p(TEVp_var = hits,
            activity = "dnAUC",
            legend = "dnAUC",
            file_name = "all_variants_dnAUC",
            height = 5, 
            asp_ratio = 1)

# heat map of variants which agree with the first uASPIre screen where
# they were identified as hits
heatmap_P1p(TEVp_var = c('CGA', 'RWG', 'LAH', 'AGM', 'ERA', 'QRS', 'DRA', 'GST', 'QSQ', 'VST', 'TSL', 'ASL', 'ATL', 'DSL'),
            activity = "dnAUC",
            legend = "dnAUC",
            file_name = "validated_in_uASPIre_dnAUC",
            height = 3.71,
            asp_ratio = 1.618)


# Comparison to in vitro data ---------------------------------------------

## prepare data ----

# define working directory
DIR <- gsub("uASPIre screen4 and 5_2022_08_11/uASPIre screen4_DMS", "", getwd())

# read in vitro data
data_vini <- read.table(paste0(DIR, "CFX in vitro assay/results/data_vini.txt"),
                          sep = "\t", header = TRUE)

# format table
data_vini <- data_vini %>%
  select(-intercept, -r_squared) %>% # remove columns not needed
  filter(!pro %in% c("buffer", "EV", "TEVp I C151A", "TEVp I VDS")) %>% # remove samples not needed
  rename("TEVp" = "pro", "TEVs_AA" = "sub") # rename columns

# rename parent variant
data_vini$TEVp[data_vini$TEVp == "TEVp I"] <- "TEVp I TSL"

# format TEVp names
data_vini <- data_vini %>% 
  mutate(TEVp_parent = substring(TEVp, 6, 6), # separate the name of TEVp into parent and position 30-32
         TEVp = substring(TEVp, 8, 10))


## to plot heat map of in vitro data ----

# write function to generate P1' heat maps
heatmap_P1p_iv <- function(
  TEVp_var, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data_vini,
                      TEVp_parent == "I",
                      TEVp %in% TEVp_var)
  
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
                             fill = vini)) +
    geom_tile() + # make heat map
    scale_fill_viridis(option="magma") +
    theme_cowplot(font_size = 12) +
    theme(axis.line = element_blank()) +
    theme(axis.title.y = element_blank()) +
    #    theme(axis.text.y = element_text(family = "Courier")) + # use a monospace font for the y-axis labels
    labs(fill = "Initial \nvelocity") +
    scale_x_discrete(name = "Amino acid in P1'", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    #    theme(axis.text.x = element_text(size = 10)) + # change size of x axis labels
    facet_grid(category ~ ., scales = "free", space = "free") # switch rows and columns
  
  # show plot
  print(p)
  # save plot to file
  save_plot(p, filename = paste0("./06_validation_of_hits/plots/heatmaps/P1p_iv_", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}

## plot iv data ----

# plot iv data
heatmap_P1p_iv(TEVp_var = pull(data_vini, TEVp),
               file_name = "all_variants",
               height = 3.71,
               asp_ratio = 1)

heatmap_P1p_iv(TEVp_var = c("TSL", "DSL", "ATL", "ASL", "GST"),
               file_name = "noGGT_resized",
               height = 3,
               asp_ratio = 1.618)


### plot iv data as bar plot ----

TEVp_var <- c("TSL", "DSL", "ATL", "ASL", "GST")

# preparing data for the plot
plot_data <- filter(data_vini,
                  TEVp_parent == "I",
                  TEVp %in% TEVp_var)

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

# as grouped bar plot
p <- ggplot(plot_data, aes(x = factor(TEVp, levels = c("GST", "TSL", "ASL", "ATL", "DSL")),
                           y = vini, fill = TEVs_AA)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  # try more space between bars
  # geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  # scale_fill_brewer(palette = "RdBu", type = "div") +  # Use the RdBu palette
  scale_fill_manual(values = c(H = "#fddbc7",
                               M = "#67a9cf",
                               S = "grey",
                               T = "#ef8a62",
                               W = "#d1e5f0")) +
  labs(x = "TEVp variant",
       y = "Initial reaction velocity in vitro",
       fill = "Amino acid in P1':") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "top")

# print plot
print(p)
# save plot
save_plot(p, filename = "./06_validation_of_hits/plots/P1p_iv_barplot_grouped.pdf",
          base_height = 3.1,
          base_asp = 1.5)


  
## plot uASPIre data again in style of iv data ----

# write function to generate P1' heat maps
heatmap_P1p_comp <- function(
  TEVp_var, TEVs_var, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data_manual,
                      TEVp_parent == "I",
                      TEVs_AA %in% TEVs_var,
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
  save_plot(p, filename = paste0("./06_validation_of_hits/plots/heatmaps/P1p_comp_", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}


# plot uASPIre data but only the TEVps and TEVs's I validated in vitro using noise-blanked values
heatmap_P1p_comp(TEVp_var = pull(data_vini, TEVp),
                 TEVs_var = pull(data_vini, TEVs_AA),
                 activity = "dnAUC",
                 legend = "dnAUC",
                 file_name = "all_variants_dnAUC",
                 height = 3.71, 
                 asp_ratio = 1)

