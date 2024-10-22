# load libraries
library(tidyverse)
library(viridis)
library(cowplot)
library(stringi)
library(readxl)

# create folders ----------------------------------------------------------

# create a folder to store the individual plots and data in
dir.create("./05_DMS")
dir.create("./05_DMS/results")
dir.create("./05_DMS/plots")
dir.create("./05_DMS/plots/heatmaps")
dir.create("./05_DMS/plots/heatmaps/orthogonal")
dir.create("./05_DMS/plots/heatmaps/promiscuous")

# prepare data ------------------------------------------------------------

# read data
data <- read.table("./04_DMS_false_positives/results/Data_DMS_fc_corrected_min20DNAlevel.txt", sep = "\t", header = TRUE)

summary_n_mut_TEVp <- data %>%
  # filter(TEVs != "ENLYFQS") %>%
  group_by(n_mut_TEVp) %>%
  summarize(avg_reads = mean(sum),
            median_reads = median(sum))

# only keep most important columns
data <- data %>% 
  select(TEVp, TEVs, n_mut_TEVp, n_mut_TEVs, min, AUC)


# put data of combinatorial TEVp mutants into a separate data frame
data_combiTEVp <- data %>% filter(n_mut_TEVp >= 2)
  
# put data of combinatorial TEVs mutants into a separate data frame
data_combiTEVs <- data %>% filter(n_mut_TEVs >= 2)

# remove all combinatorial mutants from the main data set
data <- data %>%
  filter(n_mut_TEVp <= 1,
         n_mut_TEVs <= 1)

# get single TEVs mutants only
TEVs_singlemut <- data %>% filter(n_mut_TEVs == 1)

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

# get single TEVp mutants only
TEVp_singlemut <- data %>% filter(n_mut_TEVp == 1)

# define which position in TEVp is mutated
TEVp_singlemut <- TEVp_singlemut %>% 
  mutate(TEVp_short = substring(TEVp, 2, 5)) %>%
  mutate(TEVp_rev = stri_reverse(TEVp_short)) %>%
  mutate(TEVp_pos = as.integer(stri_reverse(substring(TEVp_rev, 2, 5)))) %>%
  mutate(TEVp_AA = substring(TEVp_rev, 1, 1)) %>%
  select(-TEVp_rev, -TEVp_short)

# add the two new columns to the original data set
data <- full_join(data, TEVp_singlemut)


## Pie charts diagrams ----------------------------------------------------------

count_TEVp <- data %>%
  filter(TEVp_pos %in% seq(0, 234) | TEVp == "wt") %>%
  pull(TEVp) %>% unique() %>% length()

count_TEVp_combinatorial <- data_combiTEVp %>%
  pull(TEVp) %>% unique() %>% length()

count_TEVs <- data %>%
  filter(n_mut_TEVs <= 1) %>%
  filter(!grepl("\\*", TEVs)) %>%
  pull(TEVs) %>% unique() %>% length()

count_TEVs_combinatorial <- data_combiTEVs %>%
  pull(TEVs) %>% unique() %>% length()

count_combis <- data %>%
  filter(TEVp_pos %in% seq(0, 234) | TEVp == "wt") %>%
  filter(n_mut_TEVs <= 1) %>%
  filter(!grepl("\\*", TEVs)) %>%
  nrow()
  


cat(paste0("Percent of combinatorial TEVp mutants from all unique TEVps: ",
           round((count_TEVp_combinatorial / (count_TEVp + count_TEVp_combinatorial)*100), 2)))


# make a data frame summarizing the coverage
coverage <- data.frame(library = c("TEVp library", "TEVs library", "Combinatorial library"),
                       covered = c(count_TEVp, count_TEVs, count_combis),
                       theoretical = c(234*19+1, 7*19+1, (234*19+1) * (7*19+1)))

# calculate some numbers
coverage <- coverage %>%
  mutate(missed = theoretical - covered,
         rel_coverage = covered / theoretical)

# pivot longer for plotting
coverage_longer <- coverage %>%
  select(library, covered, missed) %>%
  pivot_longer(., cols = c(covered, missed), names_to = "category")

for (i in unique(coverage_longer$library)) {
  p <- ggplot(coverage_longer %>% filter(library == i), 
              aes(x = "", y = value, fill = category)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    labs(title = i) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank()) +
    scale_fill_manual(values = c("missed" = "#cfcfcf", "covered" = "#69498a"))
  print(p)
  # save plot
  save_plot(p, filename = paste0("./05_DMS/plots/pie_chart_", i, ".pdf"),
            base_height = 2, 
            base_asp = 1.3)
}



## dealing with noise (dnAUC, dnAUCrel)----

### correcting over all TEVs positions ----

# get information about noise
df_noise <- read.table("./03_DMS_controls/results/df_noise.txt", sep = "\t", header = TRUE)

cutoff <- df_noise %>% pull(cutoff)
mu <- df_noise %>% pull(mu)

## manipulate data acc. to noise ##

# calculate noise-blanked values (delta noise AUC, dnAUC)
# for that, set all values below the cutoff to the mean of the noise
# then blank to the mean of the noise
data$dnAUC <- ifelse(data$AUC <= cutoff, mu, data$AUC)-mu


### calculate dnAUCrel ----

# calculate noise-blanked values rel to wt TEVs (dnAUCrel)
# get values of noise-blanked wt TEVs
dnAUC_ENLYFQS <- data %>% filter(TEVs == "ENLYFQS") %>% select(TEVp, dnAUC) %>% rename("dnAUC_ENLYFQS" = "dnAUC")
# assign a very low value to cases where dnAUC of ENLYFQS is 0 as I cannot divide by 0 later
dnAUC_ENLYFQS$dnAUC_ENLYFQS[dnAUC_ENLYFQS$dnAUC_ENLYFQS == 0] <- 0.005
# join data with these values
data <- full_join(data, dnAUC_ENLYFQS, by = c("TEVp")) 
# calculate normalized blanked values
data$dnAUCrel <- data$dnAUC / data$dnAUC_ENLYFQS

# # assign the value zero where dnAUCrel is NaN
# data$dnAUCrel[is.na(data$dnAUCrel)] <- 0

# remove columns not needed
data <- data %>% select(-dnAUC_ENLYFQS)

### calculate fc_spec ----

# calculate the fold change in specificity for a given substrate compared to wt TEVp
# get values of wt TEVp
dnAUCrel_TEVpwt <- data %>% filter(TEVp == "wt") %>% select(TEVs, dnAUCrel) %>% rename("dnAUCrel_TEVpwt" = "dnAUCrel")
# find out min value besides 0
min(dnAUCrel_TEVpwt$dnAUCrel_TEVpwt[dnAUCrel_TEVpwt$dnAUCrel_TEVpwt >0], na.rm = TRUE)
# assign a very low value to cases where dnAUCrel of wt TEVp is 0 as I cannot divide by 0 later
dnAUCrel_TEVpwt$dnAUCrel_TEVpwt[dnAUCrel_TEVpwt$dnAUCrel_TEVpwt == 0] <- 0.05
# join data with these values
data <- full_join(data, dnAUCrel_TEVpwt, by = c("TEVs")) 
# calculate the fc in specificity
data$fc_spec <- data$dnAUCrel / data$dnAUCrel_TEVpwt

# # assign the value zero where fc_spec is NaN as a result of dividing by 0,
# # i.e. when dnAUCrel of wt TEVp = 0
# data$fc_spec[is.na(data$fc_spec)] <- 0

# # find out max value besides Inf values
# max(data$fc_spec[is.finite(data$fc_spec)], na.rm = TRUE)

# # assign a reasonably high value where fc_spec is Inf
# data$fc_spec[!is.finite(data$fc_spec)] <- 25

# remove columns not needed
data <- data %>% select(-dnAUCrel_TEVpwt)

## save curated data ----

# count how many unique TEVp-TEVs combinations we covered (single mutants only!)
data <- data %>%
  mutate(TEVp_TEVs = paste0(TEVp, "_", TEVs))

numbers <- length(unique(pull(data, TEVp_TEVs)))


# extract variants of the ABC DEF XYZ positions for some other checks
data_MLDEEP <- data %>% 
  filter(TEVp_pos %in% c(seq(30, 32), seq(148, 150), seq(217, 219)),
         TEVs_pos == "P1'")

# save tables to files
write.table(data_combiTEVp, file = "./05_DMS/results/combinatorial_TEVp_mutants.txt", sep = "\t", row.names=FALSE, quote = F)
write.table(data_combiTEVs, file = "./05_DMS/results/combinatorial_TEVs_mutants.txt", sep = "\t", row.names=FALSE, quote = F)
write.table(data, file = "./05_DMS/results/single_TEVp_TEVs_mutants.txt", sep = "\t", row.names=FALSE, quote = F)
write.table(data_MLDEEP, file = "./05_DMS/results/single_TEVp_TEVs_mutants_ABCDEFXYZ.txt", sep = "\t", row.names=FALSE, quote = F)

## heatmap data preparation ----

# in order to plot values of wt TEVp and wt TEVs as well in the heat maps, it is
# necessary to copy their data for each position in TEVp and TEVs

data_heatmaps <- rbind(
  # first, copy data of wildtype TEVs once for each TEVs position and append it
  # this step is necessary to avoid "holes" in heat maps
  filter(data, n_mut_TEVs == 1),
  filter(data, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P6", TEVs_AA = "E"),
  filter(data, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P5", TEVs_AA = "N"),
  filter(data, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P4", TEVs_AA = "L"),
  filter(data, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P3", TEVs_AA = "Y"),
  filter(data, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P2", TEVs_AA = "F"),
  filter(data, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P1", TEVs_AA = "Q"),
  filter(data, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P1'", TEVs_AA = "S"))


## get the wt AA of TEVp for each position
# extract the wt AA per position
TEVpwt_list <- lapply(seq(0, 235), function(pos) {
  wt_AA <- data %>%
    filter(TEVp_pos == pos) %>%
    pull(TEVp) %>%
    substr(1, 1) %>%
    unique()
  
  data.frame(TEVp = "wt",
             TEVp_pos = pos, 
             TEVp_AA = wt_AA)
})
# bind all positions together
TEVpwt <- do.call(rbind, TEVpwt_list)
# save the TEVpwt table, i.e. TEVp 0 primary structure
write.table(TEVpwt, file = "./05_DMS/results/TEVp0_seq.txt", sep = "\t", row.names=FALSE, quote = F)

# combine this data frame with actual data to generate copies of the WT activity for each position
data_WT <- left_join(TEVpwt, 
                     data_heatmaps %>% select(-TEVp_pos, -TEVp_AA), by = "TEVp")

# bind together with the rest of the data
data_heatmaps <- rbind(data_heatmaps %>% filter(TEVp != "wt"),
                       data_WT) %>%
  mutate(TEVs_pos = factor(TEVs_pos, levels=c("P1'", "P1", "P2", "P3", "P4", "P5", "P6")))

# calculate the activity normalized to the wt activity (wt TEVp with wt TEVs)
# get wt activity
wt_activity <- data %>% filter(TEVp == "wt", TEVs == "ENLYFQS") %>% pull(dnAUC)

data_heatmaps <- data_heatmaps %>%
  mutate(dnAUC_reltowtTEVpwtTEVs = dnAUC / wt_activity)

# calculate the activity normalized to the wt TEVp (wt TEVp with the respective TEVs)
# get wt TEVp activities
wtTEVp <- data %>%
  filter(TEVp == "wt") %>%
  select(TEVs, dnAUC) 

# add wtTEVp activity to heatmap data frame
data_heatmaps <- left_join(data_heatmaps, 
                  wtTEVp %>%
                    rename(dnAUC_TEVpwt = dnAUC))

# calculate activity relative to wt TEVp within each substrate
data_heatmaps <- data_heatmaps %>%
  mutate(dnAUC_reltowtTEVp = dnAUC / dnAUC_TEVpwt) %>%
  select(-dnAUC_TEVpwt)


# # save the data_heatmaps file
# write.table(data_heatmaps, file = "./05_DMS/results/data_heatmaps.txt", sep = "\t", row.names=FALSE, quote = F)

# read data
# data_heatmaps <- read.table("./05_DMS/results/data_heatmaps.txt", sep = "\t", header = TRUE, quote = "")

# save the data of TEVp 0 ("wt") only (this TEVp 0 data is more reliable than the 
# TEVp 0 data from the internal controls as it contains much more reads)
data_TEVp0 <- data_heatmaps %>% filter(TEVp == "wt")
write.table(data_TEVp0, file = "./05_DMS/results/data_TEVp0_DMS.txt", sep = "\t", row.names=FALSE, quote = F)


## get data of controls ----

# for some heat maps, we also want to plot the controls

# read data
data_controls <- read.table("./03_DMS_controls/results/data_controls.txt", sep = "\t", 
                            header = TRUE,
                            quote = "") # this is necessary to prevent R from treating the single quote in P1' as a delimiter
  
# remove unnecessary columns
data_controls <- data_controls %>%
  filter(TEVp != "0")

# in order to plot values of wt TEVs in the heat maps, it is
# necessary to copy their data for each position in TEVp and TEVs
data_controls <- rbind(
  # first, copy data of wildtype TEVs once for each TEVs position and append it
  # this step is necessary to avoid "holes" in heat maps
  filter(data_controls, n_mut_TEVs != 0),
  filter(data_controls, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P6", TEVs_AA = "E"),
  filter(data_controls, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P5", TEVs_AA = "N"),
  filter(data_controls, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P4", TEVs_AA = "L"),
  filter(data_controls, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P3", TEVs_AA = "Y"),
  filter(data_controls, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P2", TEVs_AA = "F"),
  filter(data_controls, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P1", TEVs_AA = "Q"),
  filter(data_controls, n_mut_TEVs == 0) %>% mutate(TEVs_pos = "P1'", TEVs_AA = "S"))

## get AA properties ----

# define amino acid properties for ordering them in heat maps
aa_categories <- read_excel("./data/amino_acid_properties.xlsx",
                            sheet = "classification",
                            col_names = TRUE) 

aa_categories <- pivot_longer(aa_categories, cols = names(aa_categories), names_to = "category", values_to = "AA") %>%
  filter(!is.na(AA))

aa_hydrophobicity <- read_excel("./data/amino_acid_properties.xlsx",
                                sheet = "hydrophobicity",
                                col_names = TRUE)

aa_hydrophobicity <- aa_hydrophobicity %>%
  rename("AA" = "1-letter code",
         "hydrophobicity" = "hydrophobicity index at pH7") %>%
  select(AA, hydrophobicity) %>%
  filter(!is.na(hydrophobicity))

aa_size <- read_excel("./data/amino_acid_properties.xlsx",
                      sheet = "size",
                      col_names = TRUE) %>%
  rename("AA" = "1-letter Code") %>%
  filter(!is.na(AA)) %>%
  mutate(across(c("Surface", "Volume"), as.numeric))

# make a table with all amino acid properties
aa_properties <- full_join(aa_categories, aa_hydrophobicity) %>%
  full_join(., select(aa_size, c("AA", "Surface", "Volume")))

# give the categories a shorter name
aa_properties$category[aa_properties$category == "Charged pos (basic)"] <- "basic"
aa_properties$category[aa_properties$category == "Charged neg (acidic)"] <- "acidic"
aa_properties$category[aa_properties$category == "Polar uncharged"] <- "polar"
aa_properties$category[aa_properties$category == "Special"] <- "special"
aa_properties$category[aa_properties$category == "Hydrophobic aliphatic"] <- "aliphatic"
aa_properties$category[aa_properties$category == "Hydrophobic aromatic"] <- "aromatic"

# Sanity checks -----------------------------------------------------------

## all samples ----
# Create a new column with the unique combination of TEVp and TEVs
plot_data <- data %>%
  arrange(desc(AUC)) %>%
  # filter(TEVs_AA != "*",
  #        TEVp_AA != "*") %>%
  mutate(TEVp_TEVs = interaction(TEVp, TEVs, sep = "_"),
         rank = row_number()) %>%
  mutate(category = ifelse(grepl("\\*", TEVs), "stop in TEVs", 
                           ifelse(AUC <= cutoff, "below noise", "above noise"))) 

# Create the scatter plot with AUC on the y-axis
p <- ggplot(plot_data, aes(x = rank, y = AUC, color = category)) +
  geom_point(aes(size = ifelse(TEVp_TEVs == "wt_ENLYFQ*", 5, 1)),
             alpha = 0.5) +
  geom_text(data = filter(plot_data, TEVp_TEVs == "wt_ENLYFQ*"),
            aes(label = "pos. ctrl.", size = 3), hjust = -0.2, color = "black") +
  scale_color_manual(values = c("stop in TEVs" = "#66c2a5", "below noise" = "#8da0cb", "above noise" = "#fc8d62")) +
  labs(x = "Unique TEVp-TEVs combination", y = "AUC", color = NULL) +
  theme_cowplot(font_size = 12) +
  scale_size_identity() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 0.4)) +
  # place the legend inside the plot
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

print(p)

#save plot
save_plot(p, filename = "./05_DMS/plots/scatter_AUC.png",
          base_height = 3.71, 
          base_asp = 1.2)

# count how many are above noise
n_active <- plot_data %>% filter(category == "above noise") %>% nrow()

# calculate ratio
ratio <- n_active / nrow(plot_data %>% filter(category != "stop in TEVs"))

cat("Ratio of active TEVp-TEVs combinations:", ratio*100, "%")


## truncated TEVp ----

# get wt activity
wt_activity <- data %>% filter(TEVp == "wt", TEVs == "ENLYFQS") %>% pull(dnAUC)

# plot a box plot
p <- ggplot(data %>% filter(TEVs == "ENLYFQS", TEVp_AA == "*"), 
            aes(y = dnAUC, x = factor(TEVp_pos))) +
  geom_col(fill = "#3182bd") +
  geom_point(size = 1) +
  labs(x = "Truncation (stop) at TEVp position", y = "Activity") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = wt_activity, 
             linetype = "dashed", color = "black") +
  annotate("text", x = -Inf, y = wt_activity,
           label = "TEVp 0 activity",
           hjust = -0.1, vjust = -1,
           size = 3)

print(p)

# save plot to file
save_plot(p, filename = "./05_DMS/plots/truncated_TEVp_wtsubstrate.pdf",
          base_height = 2.5, 
          base_asp = 1.5)



## MLDEEP positions ----

# define MLDEEP positions
MLDEEP_pos <- data_frame(TEVp_pos = c(30, 31, 32, 148, 149, 150, 217, 218, 219),
                         space = c('ABC', 'ABC', 'ABC', 'DEF', 'DEF', 'DEF', 'XYZ', 'XYZ', 'XYZ'))

# pick the data for the MLDEEP positions
MLDEEP_data <- data_heatmaps %>%
  filter(TEVs_pos == "P1'",
         TEVp_pos %in% MLDEEP_pos$TEVp_pos,
         TEVp_AA != "*",
         TEVs_AA != "*") %>%
  select(TEVp, TEVp_pos, AUC, fc_spec) %>%
  full_join(., MLDEEP_pos, by = "TEVp_pos")

# add all the other data to it to plot it under 'All'
combined_data <- rbind(MLDEEP_data, 
                       data %>%
                         filter(TEVs_pos == "P1'",
                                TEVp_AA != "*",
                                TEVs_AA != "*") %>%
                         select(TEVp, TEVp_pos, AUC, fc_spec) %>%
                         mutate(TEVp_pos = "All",
                                space = "complete"))

# plot a violin plot
p <- ggplot(combined_data, aes(x = factor(TEVp_pos, levels = c(30, 31, 32, 148, 149, 150, 217, 218, 219, "All")), 
                               y = AUC, 
                               fill = space)) +
  geom_violin()+
  stat_summary(
    fun = mean,
    geom = 'point',
    color = 'black') +
  labs(x = "TEVp position", y = "AUC of P1' substrates", fill = "Space") +
  theme_cowplot(font_size = 12) +
  geom_hline(yintercept = cutoff, linetype = "dashed", color = "black") +
  annotate("text", x = Inf, y = -Inf,
           label = "- - - noise cutoff",
           hjust = 1, vjust = -1,
           size = 3)

print(p)

# save plot to file
save_plot(p, filename = "./05_DMS/plots/ABC_ DEF_XYZ_activities.pdf",
          base_height = 3.71, 
          base_asp = 1.618)



# plot a violin plot of the specificity change
p <- ggplot(combined_data %>% filter(space != "DEF",
                                     TEVp_pos != "All"),
            aes(x = factor(TEVp_pos, levels = c(30, 31, 32, 217, 218, 219)), 
                y = log2(fc_spec), 
                fill = space)) +
  geom_violin()+
  stat_summary(
    fun = mean,
    geom = 'point',
    color = 'black') +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "TEVp position", y = "Fold change in P1' specificty (log2)", fill = "Cluster") +
  scale_fill_discrete(labels=c('A', 'B')) +
  theme_cowplot(font_size = 12) 

print(p)

# save plot to file
save_plot(p, filename = "./05_DMS/plots/ABC_XYZ_activities_fc_spec.pdf",
          base_height = 3.71, 
          base_asp = 1.618)


# plot a violin plot of ABC positions only
p <- ggplot(combined_data %>% filter(!space %in% c("DEF", "XYZ")),
            aes(x = factor(TEVp_pos, levels = c(30, 31, 32, "All")), 
                y = AUC, 
                fill = space)) +
  geom_violin()+
  stat_summary(
    fun = mean,
    geom = 'point',
    color = 'black') +
  labs(x = "TEVp position", y = "AUC of P1' substrates", fill = "Space") +
  theme_cowplot(font_size = 12) +
  geom_hline(yintercept = cutoff, linetype = "dashed", color = "black") +
  annotate("text", x = Inf, y = -Inf,
           label = "- - - noise cutoff",
           hjust = 1, vjust = -1,
           size = 3) +
  theme(legend.position = "none")

  
print(p)

# save plot to file
save_plot(p, filename = "./05_DMS/plots/ABC_activities.pdf",
          base_height = 3.71, 
          base_asp = 1.618)

# plot a violin plot of XZY positions only
p <- ggplot(combined_data %>% filter(space %in% c("XYZ", "complete")),
            aes(x = factor(TEVp_pos, levels = c(217, 218, 219, "All")), 
                y = AUC, 
                fill = space)) +
  geom_violin()+
  stat_summary(
    fun = mean,
    geom = 'point',
    color = 'black') +
  labs(x = "TEVp position", y = "AUC of P1' substrates", fill = "Space") +
  theme_cowplot(font_size = 12) +
  geom_hline(yintercept = cutoff, linetype = "dashed", color = "black") +
  annotate("text", x = Inf, y = -Inf,
           label = "- - - noise cutoff",
           hjust = 1, vjust = -1,
           size = 3) +
  theme(legend.position = "none")


print(p)

# save plot to file
save_plot(p, filename = "./05_DMS/plots/XYZ_activities.pdf",
          base_height = 3.71, 
          base_asp = 1.618)

## All TEVp positions ----

plot_data <- data_heatmaps %>%
  filter(TEVp_pos %in% seq(1, 234),
         # TEVp != "wt", # remove wt TEVp,
         TEVp_AA != "*",
         TEVs_AA != "*") %>%
  select(TEVp, TEVs_pos, TEVp_pos, AUC, dnAUC) %>%
  # specify categories per position to color them differently
  mutate(fraction = ifelse(.$TEVp_pos %in% seq(1, 80), 1, 
                           ifelse(.$TEVp_pos %in% seq(81, 160), 2, 3)),
         category = ifelse(.$TEVp_pos %in% c(46, 81, 151), "catalytic triad",
                           ifelse(.$TEVp_pos %in% c(30, 31, 32), "ABC",
                                  ifelse(.$TEVp_pos %in% c(217, 218, 219), "XYZ", "other"))))

# plot
p <- ggplot(plot_data, aes(x = TEVp_pos, y = AUC, group = TEVp_pos, fill = category)) +
  geom_boxplot(width = 0.8, color = "black", 
               # outlier.alpha = 0.2, outlier.size = 0.5,
               outlier.shape = NA) +
  labs(x = "TEVp position", fill = NULL) +
  theme_cowplot(font_size = 12) +
  coord_cartesian(ylim = c(0, 0.15)) +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + # make the plots turned 90 degrees
  scale_fill_manual(values = c("catalytic triad" = "#d95f02", "ABC" = "#1b9e77", "XYZ" = "#7570b3", "other" = "#f0f0f0")) +
  theme(legend.position = "top", legend.box = "horizontal", legend.justification = "center") +
  geom_hline(yintercept = cutoff, linetype = "solid", color = "#ca0020", size = 0.5) +
  geom_hline(yintercept = data %>% filter(TEVp == "wt", TEVs == "ENLYFQS") %>% pull(AUC), 
             linetype = "dashed", color = "#ca0020", size = 0.5) +
  scale_x_continuous(breaks = seq(0, 234, 10),
                     expand = c(0, 0))

p_unfaceted <- p + 
  facet_wrap(~fraction, nrow = 3, scales = "free_x") + # make the plot into two rows
  theme(strip.text = element_blank(), strip.background = element_blank()) # remove the facet labels

print(p_unfaceted)

# save plot to file
save_plot(p_unfaceted, filename = "./05_DMS/plots/all_positions_boxplot.pdf",
          base_height = 6, 
          base_asp = 1.6)


# plot activity on wt TEVs (ENLYFQS) only
plot_data <- data_heatmaps %>%
  filter(TEVp_pos %in% seq(1, 234),
         # TEVp != "wt", # remove wt TEVp,
         TEVs == "ENLYFQS",
         TEVp_AA != "*",
         TEVs_AA != "*") %>%
  select(TEVp, TEVs_pos, TEVp_pos, AUC, dnAUC) %>%
  # specify categories per position to color them differently
  mutate(fraction = ifelse(.$TEVp_pos %in% seq(1, 80), 1, 
                           ifelse(.$TEVp_pos %in% seq(81, 160), 2, 3)),
         category = ifelse(.$TEVp_pos %in% c(46, 81, 151), "catalytic triad",
                           ifelse(.$TEVp_pos %in% c(30, 31, 32), "ABC",
                                  ifelse(.$TEVp_pos %in% c(217, 218, 219), "XYZ", "other"))))

p <- ggplot(plot_data, aes(x = TEVp_pos, y = AUC, group = TEVp_pos, fill = category)) +
  geom_boxplot(width = 0.8, color = "black", outlier.shape = NA) +
  labs(x = "TEVp position", y = "AUC with wt TEVs", fill = NULL) +
  theme_cowplot(font_size = 12) +
  facet_wrap(~fraction, nrow = 3, scales = "free_x") + # make the plot into two rows
  theme(strip.text = element_blank(), strip.background = element_blank()) + # remove the facet labels
  # coord_cartesian(ylim = c(0, 0.15)) +
  coord_cartesian(ylim = c(0, 0.15)) +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + # make the plots turned 90 degrees
  scale_fill_manual(values = c("catalytic triad" = "#d95f02", "ABC" = "#1b9e77", "XYZ" = "#7570b3", "other" = "#f0f0f0")) +
  theme(legend.position = "top", legend.box = "horizontal", legend.justification = "center") +
  geom_hline(yintercept = cutoff, linetype = "solid", color = "#ca0020", size = 0.5) +
  geom_hline(yintercept = data %>% filter(TEVp == "wt", TEVs == "ENLYFQS") %>% pull(AUC), 
             linetype = "dashed", color = "#ca0020", size = 0.5) +
  scale_x_continuous(breaks = seq(0, 234, 10),
                     expand = c(0, 0))

print(p)

# save plot to file
save_plot(p, filename = "./05_DMS/plots/all_positions_boxplot_ENLYFQS_only.pdf",
          base_height = 6, 
          base_asp = 1.6)


## specificity profile with box plots ----

### plot all variants including wt TEVp as circles----

# get the data of all variants
data_p <- data_heatmaps %>%
  filter(TEVp_AA != "*",
         TEVs_AA != "*",
         TEVp_pos %in% seq(1, 234)) %>%
  # remove duplicates of wt TEVp not needed in this case
  # duplicates of wt TEVs are still needed for the facets!
  select(-TEVp_pos, -TEVp_AA) %>%
  distinct()

# add category of TEVs AA to the plot data
data_p <- full_join(data_p,
                    aa_properties %>% 
                      select(category, AA) %>%
                      rename(TEVs_AA = AA),
                    by = "TEVs_AA")

# get the data of wt TEVp
data_p_wt <- data_heatmaps %>%
  filter(TEVs_AA != "*",
         TEVp == "wt") %>%
  # remove duplicates of wt TEVp not needed in this case
  # duplicates of wt TEVs are still needed for the facets!
  select(-TEVp_pos, -TEVp_AA) %>%
  distinct()

# add category of TEVs AA to the plot data
data_p_wt <- full_join(data_p_wt,
                       aa_properties %>% 
                         select(category, AA) %>%
                         rename(TEVs_AA = AA),
                       by = "TEVs_AA")

# define the order of AA categories
category_order <- c("acidic", "basic", "aromatic", "polar", "aliphatic", "special")

# define which bars to highlight
highlight_conditions <- data.frame(
  TEVs_pos = c("P6", "P5", "P4", "P3", "P2", "P1", "P1'"),
  TEVs_AA = c("E", "N", "L", "Y", "F", "Q", "S"))

# Data for highlighting
highlight_data <- data_p %>%
  inner_join(highlight_conditions, by = c("TEVs_pos", "TEVs_AA"))


# base plot
p <- ggplot(data_p, aes(y = factor(TEVs_AA, levels = aa_properties %>%
                                     arrange(desc(hydrophobicity)) %>%
                                     arrange(desc(category)) %>%
                                     pull(AA)),
                        x = dnAUC, fill = factor(category))) +
  geom_boxplot(width = 0.8, color = "black", 
               outlier.shape = NA) +
  geom_point(data = data_p_wt, aes(y = factor(TEVs_AA, levels = aa_properties %>%
                                                arrange(desc(hydrophobicity)) %>%
                                                arrange(desc(category)) %>%
                                                pull(AA)),
                                   x = dnAUC), 
             color = "black", size = 2, shape = 1) +
  facet_grid(factor(category, levels = category_order) ~ factor(TEVs_pos, levels = rev(c("P1'", "P1", "P2", "P3", "P4", "P5", "P6"))),
             scales = "free", space = "free", switch = "y") +
  labs(y = "TEVs AA", x = "Activity") +
  theme_cowplot(font_size = 12) +
  scale_x_continuous(labels = function(x) ifelse(x == 0, "0", as.character(x))) +
  scale_fill_brewer(palette = "RdBu", direction = -1) + 
  theme(strip.placement = "outside",
        legend.position = "none")

# Data for highlighting
highlight_data <- data_p %>%
  inner_join(highlight_conditions, by = c("TEVs_pos", "TEVs_AA"))

# Overlay the highlighted geom_boxplot
p <- p + geom_boxplot(data = highlight_data, aes(y = factor(TEVs_AA),
                                                 x = dnAUC),
                      width = 0.8, color = "black", fill = NA, size = 1, outlier.shape = NA)

# show plot
print(p)
# save plot to file
save_plot(p, filename = "./05_DMS/plots/specificity_all_TEVp_mutants_horizontal_category_incl.wt.pdf",
          base_width = 9.61,
          base_asp = 3.1)



### plot wt TEVp only ----

# Data for highlighting
highlight_data <- data_p_wt %>%
  inner_join(highlight_conditions, by = c("TEVs_pos", "TEVs_AA"))

# Base ggplot
p <- ggplot(data_p_wt, aes(y = factor(TEVs_AA, levels = aa_properties %>%
                                        arrange(desc(hydrophobicity)) %>%
                                        arrange(desc(category)) %>%
                                        pull(AA)),
                           x = dnAUC, fill = factor(category))) +
  geom_col(color = "black") +
  # geom_point(size = 0.5) +
  facet_grid(factor(category, levels = category_order) ~ factor(TEVs_pos, levels = rev(c("P1'", "P1", "P2", "P3", "P4", "P5", "P6"))),
             scales = "free_y", space = "free", switch = "y") +
  labs(y = "TEVs AA", x = "Activity") +
  theme_cowplot(font_size = 12) +
  scale_x_continuous(labels = function(x) ifelse(x == 0, "0", as.character(x))) +
  scale_fill_brewer(palette = "RdBu", direction = -1) + 
  theme(strip.placement = "outside",
        legend.position = "none")

# Overlay the highlighted geom_boxplot
p <- p + geom_col(data = highlight_data, aes(y = factor(TEVs_AA),
                                             x = dnAUC),
                  color = "black", fill = NA, size = 1)

# show plot
print(p)
# save plot to file
save_plot(p, filename = "./05_DMS/plots/specificity_wt_TEVp_barplot.pdf",
          base_height = 5,
          base_width = 9.61)

### plot P92M in the same style ----

# get the data
data_p <- data_heatmaps %>%
  filter(TEVs_AA != "*",
         TEVp == "P92M") 

# add category of TEVs AA to the plot data
data_p <- full_join(data_p,
                       aa_properties %>% 
                         select(category, AA) %>%
                         rename(TEVs_AA = AA),
                       by = "TEVs_AA")

# Data for highlighting
highlight_data <- data_p %>%
  inner_join(highlight_conditions, by = c("TEVs_pos", "TEVs_AA"))

# Base ggplot
p <- ggplot(data_p, aes(y = factor(TEVs_AA, levels = aa_properties %>%
                                        arrange(desc(hydrophobicity)) %>%
                                        arrange(desc(category)) %>%
                                        pull(AA)),
                           x = dnAUC, fill = factor(category))) +
  geom_col(color = "black") +
  # geom_point(size = 0.5) +
  facet_grid(factor(category, levels = category_order) ~ factor(TEVs_pos, levels = rev(c("P1'", "P1", "P2", "P3", "P4", "P5", "P6"))),
             scales = "free_y", space = "free", switch = "y") +
  labs(y = "TEVs AA", x = "Activity") +
  theme_cowplot(font_size = 12) +
  scale_x_continuous(breaks = function(x) pretty(x, n = 3),  # Adjust `n` for the number of breaks
                     labels = function(x) ifelse(is.na(x), "n.d.", ifelse(x == 0, "0", as.character(x)))) +
  scale_fill_brewer(palette = "RdBu", direction = -1) + 
  theme(strip.placement = "outside",
        legend.position = "none")

# Overlay the highlighted geom_boxplot
p <- p + geom_col(data = highlight_data, aes(y = factor(TEVs_AA),
                                             x = dnAUC),
                  color = "black", fill = NA, size = 1)

# show plot
print(p)
# save plot to file
save_plot(p, filename = "./05_DMS/plots/specificity_P92M_TEVp_barplot.pdf",
          base_width = 9.61,
          base_asp = 3.1)

### plot all variants including wt TEVp as circles----

# get the data
data_p <- data_heatmaps %>%
  filter(TEVp_AA != "*",
         TEVs_AA != "*",
         TEVp_pos %in% seq(1, 234)) %>%
  # remove duplicates of wt TEVp not needed in this case
  # duplicates of wt TEVs are still needed for the facets!
  select(-TEVp_pos, -TEVp_AA) %>%
  distinct()

# add category of TEVs AA to the plot data
data_p <- full_join(data_p,
                    aa_properties %>% 
                      select(category, AA) %>%
                      rename(TEVs_AA = AA),
                    by = "TEVs_AA")

# define the order of AA categories
category_order <- c("acidic", "basic", "aromatic", "polar", "aliphatic", "special")

# define which bars to highlight
highlight_conditions <- data.frame(
  TEVs_pos = c("P6", "P5", "P4", "P3", "P2", "P1", "P1'"),
  TEVs_AA = c("E", "N", "L", "Y", "F", "Q", "S"))

# Data for highlighting
highlight_data <- data_p %>%
  inner_join(highlight_conditions, by = c("TEVs_pos", "TEVs_AA"))

# base plot
p <- ggplot(data_p, aes(y = factor(TEVs_AA, levels = aa_properties %>%
                                     arrange(desc(hydrophobicity)) %>%
                                     arrange(desc(category)) %>%
                                     pull(AA)),
                        x = dnAUC, fill = factor(category))) +
  geom_boxplot(width = 0.8, color = "black", 
               outlier.shape = NA) +
  geom_point(data = data_p_wt, aes(y = factor(TEVs_AA, levels = aa_properties %>%
                                                arrange(desc(hydrophobicity)) %>%
                                                arrange(desc(category)) %>%
                                                pull(AA)),
                                   x = dnAUC), 
             color = "black", size = 2, shape = 1) +
  facet_grid(factor(category, levels = category_order) ~ factor(TEVs_pos, levels = rev(c("P1'", "P1", "P2", "P3", "P4", "P5", "P6"))),
             scales = "free", space = "free", switch = "y") +
  labs(y = "TEVs AA", x = "Activity") +
  theme_cowplot(font_size = 12) +
  scale_x_continuous(labels = function(x) ifelse(x == 0, "0", as.character(x))) +
  scale_fill_brewer(palette = "RdBu", direction = -1) + 
  theme(strip.placement = "outside",
        legend.position = "none")

# Data for highlighting
highlight_data <- data_p %>%
  inner_join(highlight_conditions, by = c("TEVs_pos", "TEVs_AA"))

# Overlay the highlighted geom_boxplot
p <- p + geom_boxplot(data = highlight_data, aes(y = factor(TEVs_AA),
                                                 x = dnAUC),
                      width = 0.8, color = "black", fill = NA, size = 1, outlier.shape = NA)

# show plot
print(p)
# save plot to file
save_plot(p, filename = "./05_DMS/plots/specificity_all_TEVp_mutants_horizontal_category_incl.wt.pdf",
          base_width = 9.61,
          base_asp = 3.1)


## Get top variants for certain SOI ----

# Try to identivy top variants that are active on a certain substrate of interest (SOI)

dir.create("./05_DMS/plots/active_on_SOI")


plot_data <- data_heatmaps %>%
  filter(TEVp_pos %in% seq(1, 234),
         TEVp != "wt", # remove wt TEVp,
         TEVp_AA != "*",
         TEVs_AA != "*") %>%
  select(TEVp, TEVp_pos, TEVp_AA, TEVs_pos, TEVs_AA, dnAUC)
  
# specify fractions to positions to plot the data in 3 facets (i.e. rows)
plot_data <- plot_data %>%
  mutate(fraction = ifelse(.$TEVp_pos %in% seq(1, 80), 1, 
                           ifelse(.$TEVp_pos %in% seq(81, 160), 2, 3)))  

# add AA property
plot_data <- full_join(plot_data,
                       aa_properties %>% 
                         select(category, AA) %>%
                         rename(TEVs_AA = AA),
                       by = "TEVs_AA")

# define substrate of interest
TEVs_pos_OI <- "P3"
TEVs_AA_OI <- "C"

above3rdquartile <- plot_data %>%
  # mutate(median_dnAUC = median(dnAUC, na.rm = TRUE))
  # filter(dnAUC > median(plot_data$dnAUC, na.rm = TRUE))
  filter(dnAUC > quantile(dnAUC, probs = 0.75, na.rm = TRUE))

count <- above3rdquartile %>%
  group_by(TEVp_pos) %>%
  summarise(count = n()) %>%
  arrange(desc(count))


### as for website ----

#### identified based on 3rd quartile ----

plot_data <- data_heatmaps %>%
  filter(TEVp_pos %in% seq(1, 234),
         # TEVp != "wt", # remove wt TEVp,
         TEVp_AA != "*",
         TEVs_AA != "*") %>%
  select(TEVp, TEVp_pos, TEVp_AA, TEVs_pos, TEVs_AA, dnAUC)

#filter for substrate and most frequent positions >3rd quartile
plot_data <- plot_data %>% 
  filter(TEVs_pos == TEVs_pos_OI,
         TEVs_AA == TEVs_AA_OI,
         TEVp_pos %in% count$TEVp_pos[1:10]) # get top 10 positions

# calculate median per position
plot_data <- plot_data %>%
  group_by(TEVp_pos) %>%
  mutate(median_dnAUC = median(dnAUC, na.rm = TRUE),
         max_dnAUC = max(dnAUC)) %>%
  # add type to color wt TEVp
  mutate(type = ifelse(TEVp == "wt", "TEVp 0 (wt)", "mutant"))

p <- ggplot(plot_data, aes(x = reorder(TEVp_pos, -max_dnAUC), y = dnAUC, color = type)) +
  geom_point() +
  scale_color_manual(values = c("mutant" = "black",
                                "TEVp 0 (wt)" = "#ca0020")) +  # Apply custom colors
  labs(x = "TEVp position frequently >3rd quartile ordered by max(dnAUC)", 
       y = "Activity", 
       color = NULL) +
  theme_cowplot(font_size = 12) +
  geom_hline(yintercept = data %>% filter(TEVp == "wt", TEVs == "ENLYFQS") %>% pull(dnAUC), 
             linetype = "dashed", color = "#ca0020", linewidth = 0.5) +
  ggtitle(paste("Substrate of interest:", TEVs_AA_OI, "in", TEVs_pos_OI)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

print(p)

save_plot(p, filename = paste0("./05_DMS/plots/active_on_SOI/top_positions_3rdquartile_", TEVs_pos_OI, "_", TEVs_AA_OI, ".pdf"),
          base_height = 3.2, 
          base_asp = 1.6)

#### identified based on max activity ----

plot_data <- data_heatmaps %>%
  filter(TEVp_pos %in% seq(1, 234),
         # TEVp != "wt", # remove wt TEVp,
         TEVp_AA != "*",
         TEVs_AA != "*") %>%
  select(TEVp, TEVp_pos, TEVp_AA, TEVs_pos, TEVs_AA, dnAUC)

#filter for substrate 
plot_data <- plot_data %>% 
  filter(TEVs_pos == TEVs_pos_OI,
         TEVs_AA == TEVs_AA_OI)

#get positions with highest single activity values
top10 <- plot_data %>% 
  arrange(desc(dnAUC)) %>%
  pull(TEVp_pos) %>%
  unique() %>%
  head(10)
  
#filter for positions with top activity values
plot_data <- plot_data %>% 
  filter(TEVp_pos %in% top10)

# calculate median per position
plot_data <- plot_data %>%
  group_by(TEVp_pos) %>%
  mutate(median_dnAUC = median(dnAUC, na.rm = TRUE),
         max_dnAUC = max(dnAUC)) %>%
  # add type to color wt TEVp
  mutate(type = ifelse(TEVp == "wt", "TEVp 0 (wt)", "mutant"))

p <- ggplot(plot_data, aes(x = reorder(TEVp_pos, -max_dnAUC), y = dnAUC, color = type)) +
  geom_point() +
  scale_color_manual(values = c("mutant" = "black",
                                "TEVp 0 (wt)" = "#ca0020")) +  # Apply custom colors
  labs(x = "TEVp position with highest single activity", 
       y = "Activity", 
       color = NULL) +
  theme_cowplot(font_size = 12) +
  geom_hline(yintercept = data %>% filter(TEVp == "wt", TEVs == "ENLYFQS") %>% pull(dnAUC), 
             linetype = "dashed", color = "#ca0020", linewidth = 0.5) +
  ggtitle(paste("Substrate of interest:", TEVs_AA_OI, "in", TEVs_pos_OI)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

print(p)

save_plot(p, filename = paste0("./05_DMS/plots/active_on_SOI/top_positions_maxactivity_", TEVs_pos_OI, "_", TEVs_AA_OI, ".pdf"),
          base_height = 3.2, 
          base_asp = 1.6)


# P5 specificity ----------------------------------------------------


dir.create("./05_DMS/plots/P5analysis")


#### plot manually identified variants with potential novel P5 specificity ----

# get data
exempl <- data_heatmaps %>%
  filter(TEVs_pos == "P5",
         TEVp %in% c("wt", "G38K", "G38L", "G38T", "G38Y", "D81C", "D81K", "R101D", "R101S", "R101V", "G152T", "G152V"),
         TEVp_AA != "*",
         TEVs_AA != "*") 

# add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
# this is important to have a grey box in the heatmap for missing values
# create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations
df_combinations <- expand.grid(
  TEVp = unique(exempl$TEVp),
  TEVs_AA = unique(exempl$TEVs_AA))
# merge with data of actually measured TEVp x TEVs combinations
exempl <- exempl %>%
  full_join(., df_combinations,
            by = c("TEVp", "TEVs_AA")) %>% 
  arrange(TEVp_pos, TEVp_AA)

# make TEVp a factor for proper order in plot
exempl$TEVp <- factor(exempl$TEVp, levels = rev(unique(exempl$TEVp)))

# plotting
p <- ggplot(exempl, aes(x = TEVs_AA, 
                        y = TEVp, 
                        fill = dnAUC)) +
  geom_tile() +
  scale_fill_viridis(option="magma") +
  theme_cowplot(font_size = 12) +
  scale_x_discrete(name = "Amino acid in P5", expand = c(0, 0))+
  scale_y_discrete(name = "TEVp variant", expand = c(0, 0))+
  theme(axis.line = element_blank()) +
  labs(fill = "Activity")
  
# show plot
print(p)
# save plot to file
save_plot(p, filename = "./05_DMS/plots/P5analysis/handpicked_variants.pdf",
          base_height = 3.6, 
          base_asp = 1.3)

# Top P1 variants ---------------------------------------------------------

# get all variants tested against all 19 P1 substrates
topP1 <- data_heatmaps %>%
  filter(TEVp_AA != "*",
         TEVs != "ENLYFQS",
         TEVs_AA != "*",
         TEVs_pos == "P1")

# get top 20 variants for all 19 alternative P1 substrates
topP1_head20 <- topP1 %>%
  group_by(TEVs_AA) %>%
  arrange(desc(dnAUC)) %>%
  slice_head(n = 20) %>%
  select(TEVp, TEVs_AA) %>%
  mutate(rank = row_number())

# transform table
topP1_head20 <- topP1_head20 %>%
  pivot_wider(names_from = rank, values_from = TEVp)

# save table 
write.table(topP1_head20, file = "./05_DMS/results/top_variants_for_P1_activity.txt", sep = "\t", row.names=FALSE, quote = F)


# Heat maps   ---------------------------------------------

## b-factors wt TEVs ----

# Try to find the most important positions of TEVp for activity on the wt substrate

# For b-factors use e.g. the maximum value reached per position when tested on wt TEVs
# and set it relative to the wt activity

wt_activity <- data %>% filter(TEVp == "wt", TEVs == "ENLYFQS") %>% pull(dnAUC)

# get b_factors for ENLYFQS
bfactors_ENLYFQS <- data_heatmaps %>%
  filter(TEVp != "wt",
         TEVp_pos %in% seq(1, 234),
         TEVs_AA != "*",
         TEVp_AA != "*",
         TEVs == "ENLYFQS") %>%
  # remove copy entries of TEVp wt not used here
  select(-TEVs_pos, -TEVs_AA) %>%
  distinct() %>%
  # group by TEVp_pos and calculate some values
  group_by(TEVp_pos) %>%
  summarize(TEVs_variant = "ENLYFQS",
            tested = n(),
            above_noise = sum(dnAUC > 0),
            max = round(dnAUC[which.max(dnAUC)] / wt_activity, 3),
            max_AA = ifelse(any(dnAUC > 0), TEVp_AA[which.max(dnAUC)], "-"),
            percentile = round(quantile(dnAUC, probs = 0.9, na.rm = TRUE) / wt_activity, 3),
            avg = round(mean(dnAUC) / wt_activity, 3),
            max_log2 = log2(dnAUC[which.max(dnAUC)] / wt_activity),
            percentile_log2 = log2(
              quantile(dnAUC, probs = 0.9, na.rm = TRUE) / wt_activity),
            avg_log2 = log2(mean(dnAUC) / wt_activity),
            median_log2 = log2(median(dnAUC) / wt_activity),
            sum = sum(dnAUC))

# add wt AA to the table
bfactors_ENLYFQS <- right_join(TEVpwt %>% select(TEVp_pos, TEVp_AA) %>% rename(parent_AA = TEVp_AA),
                              bfactors_ENLYFQS)


# plot using only one specific metric, i.e. percentile_log2

# define a color palette
custom_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# define the limits for the color scale
values <- bfactors_ENLYFQS$percentile_log2
min_value <- min((values[is.finite(values)]))
max_value <- max((values[is.finite(values)]))
max_abs_value <- max(abs(min_value), abs(max_value))

# plot
p <- ggplot(bfactors_ENLYFQS, aes(x = TEVp_pos, 
                         y = TEVs_variant, 
                         fill = percentile_log2)) +
  geom_tile() +
  scale_fill_gradientn(colors = custom_palette, 
                       limits = c(-max_abs_value, max_abs_value),
                       na.value = custom_palette[1]) + # this makes sure that -Inf values resulting from 0 rel. activity (log2(0) = -Inf) are plotted with the lowest color 
  theme_cowplot(font_size = 12) +
  scale_y_discrete(name = "TEVs variant(s)", expand = c(0, 0)) +
  scale_x_continuous(name = "TEVp position", expand = c(0, 0),
                     breaks = seq(0, max(bfactors_ENLYFQS$TEVp_pos), by = 10))+
  theme(axis.line = element_blank()) +
  labs(fill = "Mutability \nscore") #+
  # labs(fill = "max. log2(fc)")
# show plot
print(p)
# save plot to file
save_plot(p, filename = "./05_DMS/plots/heatmaps/bfactors_ENLYFQS_percentile.pdf",
          base_height = 1.7, 
          base_asp = 7)


# get the positions that interact with the substrate
interactions <- read_excel("./data/TEVp_literature_variants.xlsx", 
                           sheet = "Interactions",
                           col_names = TRUE,
                           range = cell_cols(c("A", "B", "C")))

# add information of interactions with the substrate
bfactors_ENLYFQS <- full_join(bfactors_ENLYFQS,
                              interactions %>% select(TEVp_pos, Interaction))

# write table to file
write.table(bfactors_ENLYFQS, file = "./05_DMS/results/bfactors_using_wtTEVs.txt", sep = "\t", row.names=FALSE, quote = F)
# reformat table 
bfactors_reformated <- bfactors_ENLYFQS %>% arrange(percentile) %>% select(TEVp_pos, parent_AA, tested, above_noise, max, max_AA, percentile, avg, Interaction)
# write head of table to file for supplementary table
write.table(bfactors_reformated %>% head(25), file = "./05_DMS/results/bfactors_using_wtTEVs_lowest.txt", sep = "\t", row.names=FALSE, quote = F)
# write tail of table to file for supplementary table
write.table(bfactors_reformated %>% arrange(desc(percentile)) %>% head(25), file = "./05_DMS/results/bfactors_using_wtTEVs_highest.txt", sep = "\t", row.names=FALSE, quote = F)

# Calculate in how many cases I tested all 19 mutants
filtered_df <- bfactors_ENLYFQS[bfactors_ENLYFQS$tested == 19, ]
ratio <- nrow(filtered_df) / nrow(bfactors_ENLYFQS)
cat("Ratio of cases where all 19 mutants were tesed:", ratio, "\n")

### save selected b-factors for Pymol ----
pymol <- select(bfactors_ENLYFQS, percentile_log2)
# Assign a very low value to -Inf values resulting from 0 rel. activity (log2(0) = -Inf)
pymol$percentile_log2[pymol$percentile_log2 == -Inf] <- min(pymol$percentile_log2[is.finite(pymol$percentile_log2)]) *1.1
# save table
write.table(pymol, file = "./05_DMS/results/newBfactors.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


## b-factors per position ----
# This is to try to see which TEVp positions affect specificity at a certain substrate position.

# get bfactors for all non-wt TEVs variants per TEVs position
perTEVspos <- data_heatmaps %>%
  filter(TEVp != "wt",
         TEVp_pos %in% seq(1, 234),
         TEVs_AA != "*",
         TEVp_AA != "*",
         TEVs != "ENLYFQS") #%>% # this should consider all non-wt variants of TEVs

bfactors_pp <- perTEVspos %>%  
  group_by(TEVp_pos, TEVs_pos) %>%
  summarize(
    max = log2(dnAUC[which.max(dnAUC)] / wt_activity),
    percentile = log2(
      quantile(dnAUC, probs = 0.90, na.rm = TRUE) / wt_activity),
    median = log2(median(dnAUC) / wt_activity),
    mean = log2(mean(dnAUC) / wt_activity),
    sum = sum(dnAUC)) %>%
  rename(TEVs_variant = TEVs_pos)

# bind data together
# bfactors_new <- rbind(bfactors_ENLYFQS, bfactors_pp)
bfactors_new <- bfactors_pp

# define the limits for the color scale
values <- bfactors_new$percentile
min_value <- min((values[is.finite(values)]))
max_value <- max((values[is.finite(values)]))
max_abs_value <- max(abs(min_value), abs(max_value))

# plot 
p <- ggplot(bfactors_new, aes(x = TEVp_pos, 
                         y = factor(TEVs_variant, levels=c("ENLYFQS", "P1'", "P1", "P2", "P3", "P4", "P5", "P6")), 
                         fill = percentile)) +
  geom_tile() +
  scale_fill_gradientn(colors = custom_palette, 
                       # limits = c(-max_abs_value, max_abs_value),
                       na.value = custom_palette[1],
                       # na.value = "black"
                       ) + # this makes sure that -Inf values resulting from 0 rel. activity (log2(0) = -Inf) are plotted with the lowest color 
  theme_cowplot(font_size = 12) +
  scale_y_discrete(name = "TEVs variant(s)", expand = c(0, 0)) +
  scale_x_continuous(name = "TEVp position", expand = c(0, 0),
                     breaks = seq(0, max(bfactors_new$TEVp_pos), by = 10))+
  theme(axis.line = element_blank()) +
  labs(fill = "90th percentile \n of dnAUC (log2(fc))") +
  # labs(fill = "max dnAUC (log2(fc))") +
  theme(axis.text.x = element_text(size = 8)) # change size of x axis labels
# show plot
print(p)
# save plot to file
save_plot(p, filename = "./05_DMS/plots/heatmaps/bfactors_pp_90thpercentile.pdf",
          base_height = 3, 
          base_asp = 5)



## 2D DMS map: TEVp pos vs TEVs pos ----

# get wt activity
wt_activity <- data %>% filter(TEVp == "wt", TEVs == "ENLYFQS") %>% pull(dnAUC)

# define a percentile to be shown
perc <- 0.98
perc_sc <- 0.9

# preparing data for pos by pos heat maps
data_posbypos <- data_heatmaps %>%
  filter(TEVp_AA != "*",
         TEVs_AA != "*",
         TEVp != "wt", # exclude wt TEVp!
         TEVs != "ENLYFQS", # exclude wt TEVs!
         TEVp_pos %in% seq(1, 234)) %>%
  group_by(TEVp_pos, TEVs_pos) %>%
  summarize(obs = n(),
            percentile_fc_spec_rel = log2(quantile(fc_spec[!(TEVp == "wt" & TEVs == "ENLYFQS")], probs = perc_sc, na.rm = TRUE)),
            max_fc_spec_rel = log2(max(fc_spec[!(TEVp == "wt" & TEVs == "ENLYFQS")])),
            avg = mean(dnAUC),
            log2avg = log2(mean(dnAUC)),
            # when the avg relative to wt activity is calculated, it should exclude the wt activity from the avg
            avg_rel = log2(mean(dnAUC[!(TEVp == "wt" & TEVs == "ENLYFQS")]) / wt_activity),
            median = median(dnAUC),
            log2median = log2(median(dnAUC)),
            median_rel = log2(median(dnAUC[!(TEVp == "wt" & TEVs == "ENLYFQS")]) / wt_activity),
            percentile = quantile(dnAUC, probs = perc, na.rm = TRUE),
            percentile_rel = log2(quantile(dnAUC[!(TEVp == "wt" & TEVs == "ENLYFQS")], probs = perc, na.rm = TRUE) / wt_activity),
            min = min(dnAUC),
            max = max(dnAUC),
            max_rel = log2(max(dnAUC[!(TEVp == "wt" & TEVs == "ENLYFQS")]) / wt_activity),
            max_SF = max(dnAUCrel[!(TEVp == "wt" & TEVs == "ENLYFQS")]),
            percentile_SF = quantile(dnAUCrel[!(TEVp == "wt" & TEVs == "ENLYFQS")], probs = perc, na.rm = TRUE),
            percentile_SF_rel = log2(quantile(dnAUCrel[!(TEVp == "wt" & TEVs == "ENLYFQS")], probs = perc, na.rm = TRUE))) 


# write function to generate heat maps
heatmap_posbypos <- function(
  TEVspos, activity, file_name, height, asp_ratio){
  
  # prepare data for the plot
  plot_data <- data_posbypos %>%
    filter(TEVs_pos %in% TEVspos) %>%
    # # only consider cases with enough coverage
    # filter(obs >= 100) %>% 
    # give the name "values" to the column that should be plotted and is defined by "activity"
    rename(values = all_of(activity))
  
  # define a color palette for relative activities
  custom_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # define the limits for the color scale in case a relative activity is chosen
  values <- plot_data$values
  min_value <- min((values[is.finite(values)]))
  max_value <- max((values[is.finite(values)]))
  max_abs_value <- max(abs(min_value), abs(max_value))
  
  # define legend titles depending on the activity
  legend <- data_frame(
    activity = c("avg", "log2avg", "median", "log2median", "percentile", "avg_rel", "median_rel", 
                 "percentile_rel", "max", "max_rel", "percentile_SF", 
                 "percentile_SF_rel", "percentile_fc_spec_rel", "max_fc_spec_rel"),
    title = c("Avg. dnAUC", "Avg. dnAUC \n(log2)", "Median dnAUC", "Median \ndnAUC \n(log2)",
              paste0("dnAUC of the \n", perc*100, "th percentile"), 
              "log2(fc) of \navg. dnAUC", "log2(fc) of \nmedian dnAUC",
              paste0("dnAUC of the \n", perc*100, "th percentile \n(log2(fc))"),
              "Max. dnAUC", "max. dnAUC \n((log2(fc))", 
              paste0(perc*100, "th percentile of the \ndnAUC rel. to wt TEVs"),
              paste0(perc*100, "th percentile of the \ndnAUC rel. to wt TEVs \n(log2)"), 
              paste0(perc_sc*100, "th perc. of the \nspecificity change \n(log2(fc over wt TEVp)"),
              "max. specificity change \n((log2(fc))"))
  
  # plotting
  p <- ggplot(plot_data, aes(x = TEVp_pos, y = TEVs_pos, fill = values)) +
    geom_tile() + # make heat map
    theme_cowplot(font_size = 12) +
    theme(axis.line = element_blank()) +
    labs(fill = legend$title[legend$activity == activity]) +
    scale_y_discrete(name = "TEVs position", expand = c(0, 0)) +
    scale_x_continuous(name = "TEVp position", expand = c(0, 0),
                       breaks = seq(0, max(data_posbypos$TEVp_pos), by = 10)) +
    theme(legend.position = "top", legend.justification = "right",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.key.height = unit(10, "pt"),
          legend.key.width = unit(20, "pt"))
  
  # add color palette to the plot depending on whether the activity is rel or abs
  ifelse(grepl("rel", activity),
         p <- p + scale_fill_gradientn(colors = custom_palette, 
                       limits = c(-max_abs_value, max_abs_value),
                       na.value = custom_palette[1], # this makes sure that -Inf values resulting from 0 rel. activity (log2(0) = -Inf) are plotted with the lowest color
                       # na.value = "black|")
                       ),
         p <- p + scale_fill_viridis(option="magma", 
                                     na.value = ifelse(all_of(activity) %in% c("log2avg", "log2median"), magma(256)[1], magma(256)[256]))
                                     # this makes sure that Inf values resulting infinitely high values are plotted with the highest color
                                     # and -Inf values resulting in log2 transformation of cases where the avg is 0 are plotted in the lowest color
         )
  
  # show plot
  print(p)
  # save plot to file
  save_plot(p, filename = paste0("./05_DMS/plots/heatmaps/TEVp_pos_vs_TEVs_pos_", file_name, "_", activity, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}

# plot the fc of specificity compared to wt TEVp of all positions as log2
heatmap_posbypos(TEVspos = unique(data_heatmaps$TEVs_pos),
                 # TEVspos = data_heatmaps %>% filter(TEVs_pos != "P5") %>% pull(TEVs_pos) %>% unique(),
                 activity = "percentile_fc_spec_rel",
                 file_name = "complete",
                 height = 3.71, 
                 asp_ratio = 3)


# plot log2 of avg. dnAUC of all positions
heatmap_posbypos(TEVspos = unique(data_heatmaps$TEVs_pos),
                 activity = "log2avg",
                 file_name = "complete",
                 height = 2.7, 
                 asp_ratio = 3.6)


# # make plots for all different kinds of activity values to see what looks best
# for (i in c("avg", "avg_rel", "median", "log2median" ,"median_rel", "percentile", "percentile_rel", "max", "max_rel")) {
#   heatmap_posbypos(TEVspos = unique(data_heatmaps$TEVs_pos),
#                    activity = i,
#                    file_name = "complete",
#                    height = 3.71, 
#                    asp_ratio = 3)
# }


### Exemplary heat map  ----

# write a function to plot TEVp AA vs TEVs AA heatmap 
heatmap_exempl <- function(TEVspos, TEVppos){
  # get data
  exempl <- data_heatmaps %>%
    filter(TEVs_pos == TEVspos,
           TEVp_pos == TEVppos,
           TEVp_AA != "*",
           TEVs_AA != "*")
         
  # add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing values
  # create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations
  df_combinations <- expand.grid(
    TEVp_AA= unique(exempl$TEVp_AA),
    TEVs_AA = unique(exempl$TEVs_AA))
  # merge with data of actually measured TEVp x TEVs combinations
  exempl <- exempl %>%
    full_join(., df_combinations,
              by = c("TEVp_AA", "TEVs_AA"))
  
  # plotting
  p <- ggplot(exempl, aes(x = TEVp_AA, y = TEVs_AA, fill = dnAUC)) +
    geom_tile() +
    scale_fill_viridis(option="magma") +
    theme_cowplot(font_size = 12) +
    scale_y_discrete(name = paste("TEVs AA in position", TEVspos), expand = c(0, 0)) +
    scale_x_discrete(name = paste("TEVp AA in position", TEVppos), expand = c(0, 0))+
    theme(axis.line = element_blank()) +
    labs(fill = "Activity")
  
  # show plot
  print(p)
  # save plot to file
  save_plot(p, filename = paste0("./05_DMS/plots/heatmaps/exempl_", TEVspos, "_", TEVppos, ".pdf"),
            base_height = 2.9, 
            base_asp = 1.2)
}

# # check some positions that might yield specificity for P5
# for (pos in c(34, 38, 43, 47, 49, 168)) {
#   heatmap_exempl(TEVspos = "P5",
#                 TEVppos = pos)
# }

heatmap_exempl(TEVspos = "P1'",
               TEVppos = 31)

heatmap_exempl(TEVspos = "P1",
               TEVppos = 127)

heatmap_exempl(TEVspos = "P5",
               TEVppos = 46)


heatmap_exempl(TEVspos = "P5",
               TEVppos = 81)

heatmap_exempl(TEVspos = "P5",
               TEVppos = 151)

heatmap_exempl(TEVspos = "P1'",
               TEVppos = 30)

# heatmap_exempl(TEVspos = "P5",
#                TEVppos = 212)
# 
# heatmap_exempl(TEVspos = "P1",
#                TEVppos = 218)
# 
# heatmap_exempl(TEVspos = "P1",
#                TEVppos = 148)
# 
# heatmap_exempl(TEVspos = "P3",
#                TEVppos = 148)
# 
# heatmap_exempl(TEVspos = "P3",
#                TEVppos = 176)
# 
# heatmap_exempl(TEVspos = "P5",
#                TEVppos = 107)
# 
# heatmap_exempl(TEVspos = "P5",
#                TEVppos = 101)

# put several interesting positions in a grid
plot_list <- list()
count <- 1
for (pos in c(81, 151, 38, 149, 152, 165, 101, 103)) {
  # get data
  exempl <- data_heatmaps %>%
    filter(TEVs_pos == "P5",
           TEVp_pos == pos,
           TEVp_AA != "*",
           TEVs_AA != "*")
  
  # add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing value
  # create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations
  df_combinations <- expand.grid(
    TEVp_AA= unique(exempl$TEVp_AA),
    TEVs_AA = unique(exempl$TEVs_AA))
  # merge with data of actually measured TEVp x TEVs combinations
  exempl <- exempl %>%
    full_join(., df_combinations,
              by = c("TEVp_AA", "TEVs_AA"))
  
  # plotting
  p <- ggplot(exempl, aes(x = TEVp_AA, y = TEVs_AA, fill = dnAUC)) +
    geom_tile() +
    scale_fill_viridis(option="magma") +
    theme_cowplot(font_size = 12) +
    scale_y_discrete(name = "TEVs AA in position P5", expand = c(0, 0)) +
    scale_x_discrete(name = paste("TEVp AA in position", pos), expand = c(0, 0))+
    theme(axis.line = element_blank()) +
    labs(fill = "Activity")
  # save plot in a list
  plot_list[[count]] <- p
  # increase count by 1
  count <- count+1
}

# Put all plots into a grid
pg <- plot_grid(plotlist = plot_list, ncol = 2, byrow = TRUE)
pg
# save grid
save_plot(pg, filename = "./05_DMS/plots/heatmaps/interesting_positions_grid.pdf", 
          base_height = 12,
          base_asp = 0.8)

# put several interesting positions in a grid
plot_list <- list()
count <- 1
for (pos in c(34, 38, 46, 48, 81, 101, 103, 107, 149, 151, 152, 165, 167, 169, 178)) {
  # get data
  exempl <- data_heatmaps %>%
    filter(TEVs_pos == "P5",
           TEVp_pos == pos,
           TEVp_AA != "*",
           TEVs_AA != "*")
  
  # add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing values
  # create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations
  df_combinations <- expand.grid(
    TEVp_AA= unique(exempl$TEVp_AA),
    TEVs_AA = unique(exempl$TEVs_AA))
  # merge with data of actually measured TEVp x TEVs combinations
  exempl <- exempl %>%
    full_join(., df_combinations,
              by = c("TEVp_AA", "TEVs_AA"))
  
  # plotting
  p <- ggplot(exempl, aes(x = TEVp_AA, y = TEVs_AA, fill = dnAUC)) +
    geom_tile() +
    scale_fill_viridis(option="magma", 
                       # na.value = magma(256)[1]
                       # limits = c(0, 0.1)
                       ) +
    theme_cowplot(font_size = 10) +
    scale_y_discrete(name = "TEVs AA in position P5", expand = c(0, 0)) +
    scale_x_discrete(name = paste("TEVp AA in position", pos), expand = c(0, 0))+
    theme(axis.line = element_blank()) +
    # guides(fill = 'none')
    labs(fill = "Activity")
  # save plot in a list
  plot_list[[count]] <- p
  # increase count by 1
  count <- count+1
}

# Put all plots into a grid
pg <- plot_grid(plotlist = plot_list, ncol = 3, byrow = TRUE, labels = "auto")
pg
# save grid
save_plot(pg, filename = "./05_DMS/plots/heatmaps/interesting_positions2_grid_withlegend_smallertext.pdf", 
          base_height = 13,
          base_asp = 0.75)

# Ratios that work well:
# base_height = 12,
# base_asp = 0.8)
# 
# base_height = 14,
# base_asp = 0.65)

## Exempl. TEVs pos x TEVs AA ----

# write a function to generate the plots
heatmap_TEVsposbyTEVsAA <- function(activity, TEVspos, TEVppos){
  # get data
  exempl <- data_heatmaps %>%
    filter(TEVs_pos %in% TEVspos,
           TEVp_pos == TEVppos,
           TEVp_AA != "*",
           TEVs_AA != "*") %>%
    group_by(TEVs_pos, TEVs_AA) %>%
    summarize(avg = mean(dnAUC),
              max = dnAUC[which.max(dnAUC)],
              AAmax = TEVs_AA[which.max(dnAUC)]) %>%
    # give the name "values" to the column that should be plotted and is defined by "activity"
    rename(values = all_of(activity))
  
  # define legend titles depending on the activity
  legend <- data_frame(
    activity = c("avg", "max"),
    title = c("Avg. dnAUC", "Max. dnAUC"))
  
  # plotting
  p <- ggplot(exempl, aes(x = TEVs_AA, y = TEVs_pos, fill = values)) +
    geom_tile() + # make heat map
    theme_cowplot(font_size = 12) +
    scale_fill_viridis(option="magma") +
    theme(axis.line = element_blank()) +
    labs(fill = legend$title[legend$activity == activity]) +
    scale_y_discrete(name = "TEVs position", expand = c(0, 0)) +
    scale_x_discrete(name = paste0("TEVs AA", " (TEVp pos = ", TEVppos, ")"), expand = c(0, 0))
  
  # show plot
  print(p)
  
  # save plot to file
  save_plot(p, filename = paste0("./05_DMS/plots/heatmaps/exempl_TEVsposbyTEVsAA_", TEVppos, "_", activity,".pdf"),
            base_height = 3.71, 
            base_asp = 1.2)
}

# plot examples
heatmap_TEVsposbyTEVsAA(activity = "avg",
                        TEVspos = unique(data_heatmaps$TEVs_pos),
                        TEVppos = 46)


## Exempl. TEVs pos x TEVp AA ----

# write a function to generate the plots
heatmap_TEVsposbyTEVpAA <- function(activity, TEVspos, TEVppos){
  # get data
  exempl <- data_heatmaps %>%
    filter(TEVs_pos %in% TEVspos,
           TEVp_pos == TEVppos,
           TEVp_AA != "*",
           TEVs_AA != "*") %>%
    group_by(TEVs_pos, TEVp_AA) %>%
    summarize(avg = mean(dnAUC),
              max = dnAUC[which.max(dnAUC)],
              AAmax = TEVs_AA[which.max(dnAUC)]) %>%
    # give the name "values" to the column that should be plotted and is defined by "activity"
    rename(values = all_of(activity))
  
  # define legend titles depending on the activity
  legend <- data_frame(
    activity = c("avg", "max"),
    title = c("Avg. dnAUC", "Max. dnAUC"))
  
  # plotting
  p <- ggplot(exempl, aes(x = TEVp_AA, y = TEVs_pos, fill = values)) +
    geom_tile() + # make heat map
    theme_cowplot(font_size = 12) +
    scale_fill_viridis(option="magma") +
    theme(axis.line = element_blank()) +
    labs(fill = legend$title[legend$activity == activity]) +
    scale_y_discrete(name = "TEVs position", expand = c(0, 0)) +
    scale_x_discrete(name = paste("TEVp AA in position", TEVppos), expand = c(0, 0))
  
  # show plot
  print(p)
  
  # save plot to file
  save_plot(p, filename = paste0("./05_DMS/plots/heatmaps/exempl_TEVsposbyTEVpAA_", TEVppos, "_", activity,".pdf"),
            base_height = 3.71, 
            base_asp = 1.2)
}

# plot examples
heatmap_TEVsposbyTEVpAA(activity = "avg",
                        TEVspos = unique(data_heatmaps$TEVs_pos),
                        TEVppos = 46)


## TEVp pos vs TEVp AA for wt TEVs only----

data_p <- data_heatmaps %>%
  filter(TEVs == "ENLYFQS",
         TEVp_pos %in% seq(1, 234),
         TEVs_AA != "*",
         TEVp_AA != "*") %>%
  # remove copy entries of TEVp wt not used here
  select(-TEVs_pos, -TEVs_AA) %>%
  distinct()

# add combinations that were not covered in the screen and set the value to NA
# this is important to have a grey box in the heatmap for missing values
# create mock data frame with all unique combinations
df_combinations <- expand.grid(
  TEVp_AA = unique(data_p$TEVp_AA),
  TEVp_pos = unique(data_p$TEVp_pos))
# merge with data of actually measured combinations
data_p <- data_p %>%
  full_join(., df_combinations,
            by = c("TEVp_AA", "TEVp_pos"))

# merge with AA properties
data_p <- full_join(data_p, 
                    aa_properties %>% select(AA, category) %>% rename("TEVp_AA" = "AA"), 
                    by = "TEVp_AA")

# define the order of AA categories
category_order <- c("acidic", "basic", "aromatic", "polar", "aliphatic", "special")

# Create a named vector for the labels
x_labels <- setNames(as.character(filter(TEVpwt, TEVp_pos %in% 1:234) %>% 
                                    pull(TEVp_AA)),
                     1:234)


# calculate the log2 of dnAUC
data_p$log2 <- log2(data_p$dnAUC)
# manually assign the lowest value to all -Inf cases resulting from a dnAUC of zero
# This is necessary to plot those tiles with the lowest color
data_p$log2[data_p$log2 == -Inf] <- min(data_p$log2[is.finite(data_p$log2)], na.rm = TRUE)


#plot again using log2
p <- ggplot(data_p, aes(x = TEVp_pos, 
                        y = factor(TEVp_AA, levels = aa_properties %>% arrange(desc(hydrophobicity)) %>% pull(AA)),
                        fill = log2)) +
  geom_tile() +
  facet_grid(factor(category, levels = category_order) ~ ., scales = "free", space = "free") +
  scale_fill_viridis(option="magma") +
  theme_cowplot(font_size = 12) +
  scale_y_discrete(name = "TEVp AA", expand = c(0, 0)) +
  scale_x_continuous(name = "TEVp position", expand = c(0, 0),
                     breaks = seq(0, 234, by = 10))+
  theme(axis.line = element_blank()) +
  labs(fill = "Activity \n(log2)") +
  # theme(panel.spacing.y = unit(0.05, "lines"))  # change distance between facets
  theme(panel.spacing.y = unit(0.15, "lines")) + # change distance between facets
  theme(legend.position = "top", legend.justification = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.height = unit(10, "pt"),
        legend.key.width = unit(20, "pt"))

  
# show plot
print(p)
# save plot to file
save_plot(p, filename = "./05_DMS/plots/heatmaps/TEVp_pos_vs_TEVp_AA_wtTEVs_facet_hydrophobicity_log2.pdf",
          base_height = 3.9, 
          base_asp = 2.5)


###plot using log2norm with blue-white-red scale

# calculate the log2 of dnAUC_reltowtTEVpwtTEVs
data_p$log2norm <- log2(data_p$dnAUC_reltowtTEVpwtTEVs)
# manually assign the lowest value to all -Inf cases resulting from a dnAUC of zero
# This is necessary to plot those tiles with the lowest color
data_p$log2norm[data_p$log2norm == -Inf] <- min(data_p$log2norm[is.finite(data_p$log2norm)], na.rm = TRUE)

# define a color palette
custom_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# define the limits for the color scale
values <- data_p$log2norm
min_value <- min((values[is.finite(values)]))
max_value <- max((values[is.finite(values)]))
max_abs_value <- max(abs(min_value), abs(max_value))

p <- ggplot(data_p, aes(x = TEVp_pos, 
                        y = factor(TEVp_AA, levels = aa_properties %>% arrange(desc(hydrophobicity)) %>% pull(AA)),
                        fill = log2norm)) +
  geom_tile() +
  facet_grid(factor(category, levels = category_order) ~ ., scales = "free", space = "free") +
  # scale_fill_viridis(option="magma") +
  scale_fill_gradientn(colors = custom_palette, 
                       limits = c(-max_abs_value, max_abs_value)) +
  theme_cowplot(font_size = 12) +
  scale_y_discrete(name = "TEVp AA", expand = c(0, 0)) +
  scale_x_continuous(name = "TEVp position", expand = c(0, 0),
                     breaks = seq(0, 234, by = 10))+
  # scale_x_discrete(labels = x_labels) +
  theme(axis.line = element_blank()) +
  labs(fill = "Activity relative to \nwt TEVp (log2)") +
  # theme(panel.spacing.y = unit(0.05, "lines"))  # change distance between facets
  theme(panel.spacing.y = unit(0.15, "lines")) + # change distance between facets
  theme(legend.position = "top", legend.justification = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.height = unit(10, "pt"),
        legend.key.width = unit(20, "pt"))


# show plot
print(p)
# save plot to file
save_plot(p, filename = "./05_DMS/plots/heatmaps/TEVp_pos_vs_TEVp_AA_wtTEVs_facet_hydrophobicity_log2_normalized_blue-red.pdf",
          base_height = 3.9, 
          base_asp = 2.5)

#plot again over two rows
library(gridExtra)

part1 <- ggplot(data_p %>% filter(TEVp_pos %in% 1:117), aes(x = TEVp_pos, 
                                                            y = factor(TEVp_AA, levels = aa_properties %>% arrange(desc(hydrophobicity)) %>% pull(AA)),
                                                            fill = log2norm)) +
  geom_tile() +
  facet_grid(factor(category, levels = category_order) ~ ., scales = "free", space = "free") +
  # scale_fill_viridis(option="magma") +
  scale_fill_gradientn(colors = custom_palette, 
                       limits = c(-max_abs_value, max_abs_value)) +
  theme_cowplot(font_size = 12) +
  scale_y_discrete(name = "TEVp AA", expand = c(0, 0)) +
  # scale_x_continuous(name = "TEVp position", expand = c(0, 0),
  #                    breaks = seq(0, 234, by = 10))+
  scale_x_continuous(name = "Parent residue", expand = c(0, 0),
                     breaks = as.numeric(names(x_labels)), labels = x_labels) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.line = element_blank()) +
  labs(fill = "Activity relative to \nwt TEVp (log2)") +
  theme(panel.spacing.y = unit(0.15, "lines")) + # change distance between facets
  theme(legend.position = "top", legend.justification = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.height = unit(10, "pt"),
        legend.key.width = unit(20, "pt"))

part2 <- ggplot(data_p %>% filter(TEVp_pos %in% 118:234), aes(x = TEVp_pos, 
                                                              y = factor(TEVp_AA, levels = aa_properties %>% arrange(desc(hydrophobicity)) %>% pull(AA)),
                                                              fill = log2norm)) +
  geom_tile() +
  facet_grid(factor(category, levels = category_order) ~ ., scales = "free", space = "free") +
  # scale_fill_viridis(option="magma") +
  scale_fill_gradientn(colors = custom_palette, 
                       limits = c(-max_abs_value, max_abs_value)) +
  theme_cowplot(font_size = 12) +
  scale_y_discrete(name = "TEVp AA", expand = c(0, 0)) +
  # scale_x_continuous(name = "TEVp position", expand = c(0, 0),
  #                    breaks = seq(0, 234, by = 10))+
  scale_x_continuous(name = "Parent residue", expand = c(0, 0),
                     breaks = as.numeric(names(x_labels)), labels = x_labels) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.line = element_blank()) +
  labs(fill = "Activity relative to \nwt TEVp (log2)") +
  theme(panel.spacing.y = unit(0.15, "lines")) + # change distance between facets
  theme(legend.position = "top", legend.justification = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.height = unit(10, "pt"),
        legend.key.width = unit(20, "pt"))

p <- grid.arrange(
  part1,
  part2,
  nrow = 2)

# show plot
print(p)
# save plot to file
save_plot(p, filename = "./05_DMS/plots/heatmaps/TEVp_pos_vs_TEVp_AA_wtTEVs_facet_hydrophobicity_log2_normalized_split_blue-red.pdf",
          base_width = 9.75, 
          base_height = 7.9)

# calculate the coverage

# Count the number of combinations covered, i.e. where we have a dnAUC value
covered <- sum(!is.na(data_p$dnAUC))
missed <- sum(is.na(data_p$dnAUC))

coverage <- covered / (covered + missed)

cat("TEVp variants covered:", covered)
cat("TEVp variants missed:", missed)
cat("Coverage:", coverage)



# 3_identification of interesting variants ----------------------------------

## orthogonal variants ----

# to identify orthogonal variants divide the first maximum by the second maximum
# within each TEVp
orthogonality <- data %>%
  filter(TEVp_pos %in% seq(1, 234),
         TEVp_AA != "*", 
         TEVs_AA != "*") %>%
  group_by(TEVp, TEVs_pos) %>%
  arrange(desc(dnAUC)) %>%
  summarize(max1_AUC = max(AUC),
            max2_AUC = nth(AUC, 2),
            orthog = max1_AUC / max2_AUC,
            top_TEVs_AA = TEVs_AA[which.max(AUC)],
            substrates = n()) %>%
  arrange(desc(orthog))

# save data table
write.table(orthogonality, file = "./05_DMS/results/orthogonality.txt", sep = "\t", row.names=FALSE, quote = F)

## promiscuous variants ----

# install and load the entropy package 
library(entropy)

# to identify promiscuous variants calculate the entropy of AUC within each TEVp
# and multiply it by the mean activity of that TEVp
promiscuity <- data %>%
  filter(TEVp_pos %in% seq(1, 234),
         TEVp_AA != "*", 
         TEVs_AA != "*") %>%
  group_by(TEVp, TEVp_pos, TEVp_AA, TEVs_pos) %>%
  summarize(promisc = entropy(AUC) * mean(AUC),
            substrates = n()) %>%
  arrange(desc(promisc))


# save data table
write.table(promiscuity, file = "./05_DMS/results/promiscuity.txt", sep = "\t", row.names=FALSE, quote = F)

## to plot interesting variants ----

# write function to generate heat maps of promiscuous hits
heatmap_promiscuous <- function(
  TEVspos, TEVp_var, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data_heatmaps,
                      TEVs_pos == TEVspos,
                      TEVp %in% TEVp_var | TEVp == "wt",
                      !TEVp %in% c("E222L", "F116C"), #remove 2 variants that show suspiciously high P1 activity
                      TEVs_AA != "*") %>%
    # remove copy entries of TEVp wt not used here
    select(-TEVp_pos, -TEVp_AA) %>%
    distinct()
  
  # add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing values
  # create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations - 3381 lines
  df_combinations <- expand.grid(
    TEVp = unique(plot_data$TEVp),
    TEVs_AA = unique(aa_properties$AA))
  # merge with data of actually measured TEVp x TEVs combinations
  plot_data <- plot_data %>%
    full_join(., df_combinations,
              by = c("TEVp", "TEVs_AA")) %>%
    # give the name "values" to the column that should be plotted and is defined by "activity"
    rename(values = all_of(activity))
  
  # manually define the order of TEVp variants
  TEVp_levels <- plot_data %>%
    filter(TEVp != "wt") %>%
    group_by(TEVp) %>% 
    summarize(activity_mean = mean(na.omit(values))) %>% 
    arrange(desc(activity_mean)) %>% 
    pull(TEVp)
  
  # wt TEVp should be in the first position of the TEVp order
  TEVp_levels <- c("wt", TEVp_levels)
  
  # plotting
  p <- ggplot(plot_data, aes(x = factor(TEVp, levels = TEVp_levels), 
                             y = factor(TEVs_AA, levels = rev(sort(unique(TEVs_AA)))),
                             fill = values)) +
    geom_tile() + # make heat map
    scale_fill_viridis(option="magma") +
    theme_cowplot(font_size = 12) +
    theme(axis.line = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    xlab("") +
    labs(fill = legend) +
    scale_y_discrete(name = paste("Amino acid in", TEVspos), expand = c(0, 0))
  
  # show plot
  print(p)
  
  # save plot to file
  save_plot(p, filename = paste0("./05_DMS/plots/heatmaps/promiscuous/", TEVspos, "_", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}

# write function to generate heat maps of orthogonal hits
heatmap_orthogonal <- function(
  TEVspos, TEVp_var, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- filter(data_heatmaps,
                      TEVs_pos == TEVspos,
                      TEVp %in% TEVp_var | TEVp == "wt",
                      TEVs_AA != "*") %>%
    # remove copy entries of TEVp wt not used here
    select(-TEVp_pos, -TEVp_AA) %>%
    distinct()
  
  # add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing values
  # create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations - 3381 lines
  df_combinations <- expand.grid(
    TEVp = unique(plot_data$TEVp),
    TEVs_AA = unique(aa_properties$AA))
  # merge with data of actually measured TEVp x TEVs combinations
  plot_data <- plot_data %>%
    full_join(., df_combinations,
              by = c("TEVp", "TEVs_AA")) %>%
    # give the name "values" to the column that should be plotted and is defined by "activity"
    rename(values = all_of(activity))
  
  # manually define the order of TEVp variants, order acc. to activity of top TEVs_AA
  TEVp_levels <- plot_data %>%
    filter(TEVp != "wt") %>%
    arrange(desc(values)) %>% 
    pull(TEVp) %>% unique()
  
  # wt TEVp should be in the first position of the TEVp order
  TEVp_levels <- c("wt", TEVp_levels)
  
  # plotting
  p <- ggplot(plot_data, aes(x = factor(TEVp, levels = TEVp_levels), 
                             y = factor(TEVs_AA, levels = rev(sort(unique(TEVs_AA)))),
                             fill = values)) +
    geom_tile() + # make heat map
    scale_fill_viridis(option="magma", limits = c(0, 
                                                  data_heatmaps %>% 
                                                    filter(TEVs_AA != "*") %>% 
                                                    pull(all_of(activity)) %>% 
                                                    max())) +
    theme_cowplot(font_size = 12) +
    theme(axis.line = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    xlab("") +
    labs(fill = legend) +
    scale_y_discrete(name = paste("Amino acid in", TEVspos), expand = c(0, 0))
  
  # show plot
  print(p)
  
  # save plot to file
  save_plot(p, filename = paste0("./05_DMS/plots/heatmaps/orthogonal/", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}

# write function to generate a full specificity profile
heatmap_full_spec <- function(
    TEVp_var, activity, legend, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- data_heatmaps %>%
    # add TEVp controls
    full_join(., data_controls) %>%
    # filter for TEVp variants
    filter(TEVs_AA != "*", 
           TEVp %in% TEVp_var) %>%
    # remove copy entries of wt TEVp not needed here
    select(-TEVp_AA, -TEVp_pos) %>%
    distinct()
  
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
    facet_wrap(~TEVp) +
    scale_fill_viridis(option="magma") +
    theme_cowplot(font_size = 12) +
    theme(axis.line = element_blank()) +
    scale_y_discrete(name = "TEVs position", expand = c(0, 0)) +
    scale_x_discrete(name = "Amino acid", expand = c(0, 0)) +
    labs(fill = legend)
  
  # show plot
  print(p)
  # save plot to file
  save_plot(p, filename = paste0("./05_DMS/plots/heatmaps/", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}


## plotting of interesting variants ----

# heat map of top promiscuous variants for each TEVs position
for (i in c("P6",  "P5",  "P4",  "P3",  "P2",  "P1", "P1'")) {
  heatmap_promiscuous(TEVspos = i,
                      TEVp_var = promiscuity %>%
                        filter(TEVs_pos == i, 
                               substrates >= 15) %>%
                        arrange(desc(promisc)) %>%
                        head(50) %>% 
                        pull(TEVp),
                      activity = "dnAUC",
                      legend = "Activity",
                      file_name = "15substrates_dnAUC", 
                      height = 3.4, 
                      asp_ratio = 2.8)
}



# heat map of top orthogonal variants (irrespective which AA they are orthogonal for)
for (i in c("P6",  "P5",  "P4",  "P3",  "P2",  "P1", "P1'")) {
    heatmap_orthogonal(TEVspos = i,
                     TEVp_var = orthogonality %>%
                       filter(TEVs_pos == i,
                              substrates >= 15) %>%
                       arrange(desc(orthog)) %>%
                       head(50) %>% 
                       pull(TEVp),
                     activity = "dnAUC",
                     legend = "dnAUC",
                     file_name = paste0(i, "_15substrates_dnAUC_best"), 
                     height = 4.15, 
                     asp_ratio = 2.5)
}

## TEVp pos vs TEVs AA per TEVs pos ----

### Avg. activity ----

# initialize an empty list to store the plots in
plot_list <- list()

# define a threshold for how many TEVp mutants should have been tested per TEVs
threshold <- 5

for (TEVspos in c("P6",  "P5",  "P4",  "P3",  "P2",  "P1", "P1'")) {
  
  data_posbyTEVsAA <- data_heatmaps %>%
    filter(TEVs_pos == TEVspos,
           TEVp_AA != "*",
           TEVs_AA != "*",
           !TEVp_pos %in% c(0, 235)) %>%
    group_by(TEVp_pos, TEVs_AA) %>%
    summarize(avg_dnAUC = mean(dnAUC),
              min = min(dnAUC),
              max = max(dnAUC),
              AAmax = TEVp_AA[which.max(dnAUC)],
              obs = n(),
              obs_rel = n()/20) %>%
    filter(obs >= threshold)
  
  # add the wt AA to the data frame
  data_posbyTEVsAA <- left_join(data_posbyTEVsAA, TEVpwt, by = "TEVp_pos") %>%
    select(-TEVp) %>%
    rename(AAwt = TEVp_AA)
  
  # add combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing values
  # create mock data frame with all unique combinations
  df_combinations <- expand.grid(
    TEVs_AA = unique(data_posbyTEVsAA$TEVs_AA),
    TEVp_pos = unique(data_posbyTEVsAA$TEVp_pos))
  # merge with data of actually measured combinations
  data_p <- data_posbyTEVsAA %>%
    full_join(., df_combinations,
              by = c("TEVs_AA", "TEVp_pos"))
  
  # plotting
  p <- ggplot(data_p, aes(x = TEVp_pos, y = TEVs_AA, fill = avg_dnAUC)) +
    geom_tile() +
    scale_fill_viridis(option="magma") +
    theme_cowplot(font_size = 12) +
    scale_y_discrete(name = "TEVs AA", expand = c(0, 0)) +
    scale_x_continuous(name = "TEVp position", expand = c(0, 0),
                       breaks = seq(0, 234, by = 10))+
    theme(axis.line = element_blank()) +
    # theme(axis.title.y = element_blank()) +
    labs(fill = paste("Avg. \nactivity \nin", TEVspos)) +
    theme(axis.text.x = element_text(size = 8)) # change size of x axis labels
  
  # show plot
  print(p)
  # save plot to file
  save_plot(p, filename = paste0("./05_DMS/plots/heatmaps/TEVp_pos_vs_TEVs_AA_in_", TEVspos, ".pdf"),
            base_height = 3.71, 
            base_asp = 3)
  
  # save plot in a list
  plot_list[[TEVspos]] <- p
}

# Put all plots into a grid
pg <- plot_grid(plotlist = plot_list, ncol = 1, byrow = FALSE)

# save grid
save_plot(pg, filename = "./05_DMS/plots/heatmaps/TEVp_pos_vs_TEVs_AA_per_TEVspos_grid.pdf", 
          base_height = 19,
          base_asp = 0.6)


## Looking for desired activity ----

TEVspos <- "P1'"

data_posbyTEVsAA <- data_heatmaps %>%
  filter(TEVs_pos == TEVspos,
         TEVp_AA != "*",
         TEVs_AA != "*",
         TEVp_pos %in% seq(1, 234)) %>%
  group_by(TEVp_pos, TEVs_AA) %>%
  summarize(avg_dnAUC = mean(dnAUC),
            min = min(dnAUC),
            max = max(dnAUC),
            AAmax = TEVp_AA[which.max(dnAUC)],
            obs = n(),
            obs_rel = n()/20)

# add the wt AA to the data frame
data_posbyTEVsAA <- left_join(data_posbyTEVsAA, TEVpwt, by = "TEVp_pos") %>%
  select(-TEVp) %>%
  rename(AAwt = TEVp_AA)


max(data_posbyTEVsAA$obs)

# define an exemplary substrate of interest
TEVspos <- "P1'"
TEVsAA <- "N"

# get the data
test <- data_heatmaps %>%
  filter(TEVs_pos == TEVspos,
         TEVs_AA == TEVsAA,
         TEVp_pos %in% seq(1, 234),
         TEVs_AA != "*",
         TEVp_AA != "*")

# add combinations that were not covered in the screen and set the value to NA
# this is important to have a grey box in the heatmap for missing values
# create mock data frame with all unique combinations
df_combinations <- expand.grid(
  TEVp_AA = unique(test$TEVp_AA),
  TEVp_pos = unique(test$TEVp_pos))
# merge with data of actually measured combinations
data_p <- test %>%
  full_join(., df_combinations,
            by = c("TEVp_AA", "TEVp_pos"))

# plotting
p <- ggplot(data_p, aes(x = TEVp_pos, y = TEVp_AA, fill = dnAUC)) +
  geom_tile() +
  scale_fill_viridis(option="magma") +
  theme_cowplot(font_size = 12) +
  scale_y_discrete(name = "TEVp AA", expand = c(0, 0)) +
  scale_x_continuous(name = "TEVp position", expand = c(0, 0),
                     breaks = seq(0, 234, by = 10))+
  theme(axis.line = element_blank()) +
  labs(fill = paste0("dnAUC using\n", TEVsAA, " in ", TEVspos)) +
  # theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 8)) # change size of x axis labels

# show plot
print(p)
# save plot to file
save_plot(p, filename = paste0("./05_DMS/plots/heatmaps/exempl_TEVp_pos_vs_TEVp_AA_", TEVspos, TEVsAA, ".pdf"),
          base_height = 3.71, 
          base_asp = 3)


# Literature variants -----------------------------------------------------

# What follows are comparisons of our data with TEVp variants reported in the 
# literature to have altered substrate specificity.

## get literature variants ----
# get variants from literature
literature_variants <- read_excel("./data/TEVp_literature_variants.xlsx", 
                                  sheet = "Sheet1",
                                  col_names = TRUE) %>% 
  head(27) %>%
  select(-Library_generation, -Note, -Assumed_MOA)

# bring into longer format
literature_variants_long <- pivot_longer(literature_variants,
                                         cols = paste("Mutation", seq(1:20)),
                                         names_to = 'mutation_no',
                                         values_to = "mutation") %>%
  select(-mutation_no) %>%
  filter(!is.na(mutation))

# add a column specifying the mutated TEVp position
literature_variants_long <- literature_variants_long %>%
  mutate(TEVp_pos = as.numeric(str_sub(mutation, start = 2, end = -2))) 


## to plot literature variants ----

# generate a folder to store the plots in
dir.create("./05_DMS/plots/literature_variants")

# write function to generate plots investigating the individual mutations of 
# variants from literature
check_plot <- function(
  variant, TEVspos, TEVsAA, activity, height, asp_ratio){
  # get the variant's mutations and mutated TEVp pos
  variant_info <- literature_variants_long %>%
    filter(Variant_name == variant) %>%
    mutate(variantAA = str_sub(mutation, -1))
  
  # get relevant data
  exempl <- data_heatmaps %>%
    filter(TEVp_pos %in% variant_info$TEVp_pos,
           TEVs_pos == TEVspos,
           TEVs_AA == TEVsAA,
           TEVp_AA != "*") %>%
    # add information which TEVp_AA yields the highest dnAUC
    group_by(TEVp_pos) %>%
    mutate(AAmax = TEVp_AA[which.max(dnAUC)]) %>%
    # make a separate column stating the wt AA per position
    ungroup %>%
    mutate(AAwt = substr(.$TEVp, 1, 1))
  
  # add the column of the variant's AA to this data
  exempl <- left_join(exempl,
                      variant_info %>% select(TEVp_pos, variantAA),
                      by = "TEVp_pos")
  
  # give the inactive cases a very small activity, so they can be 
  # distinguished from cases that were not tested
  exempl$dnAUC[exempl$dnAUC == 0] <- 0.001
  
  # add all TEVp AAs, even if they were not covered
  df_combinations <- expand.grid(
    TEVp_AA = unique(aa_properties$AA),
    TEVp_pos = unique(exempl$TEVp_pos))
  # merge with data of actually measured combinations
  data_p <- full_join(exempl,
                      df_combinations,
                      by = c("TEVp_AA", "TEVp_pos")) %>%
    # give the name "values" to the column that should be plotted and is defined by "activity"
    rename(values = all_of(activity))
  
  
  # make a data frame with labelling informations
  variant_info <- full_join(variant_info,
                            exempl %>%
                              filter(TEVp %in% variant_info$mutation) %>%
                              select(TEVp, dnAUC, AAmax, AAwt) %>%
                              rename(mutation = TEVp)) %>%
    mutate(TEVp_AA = str_sub(mutation, start = -1, end = -1))
  
  # define the max and middle activity value to display
  max_value <- round(max(data_p$values, na.rm=T), 2)
  middle_value <- max_value/2
  
  # plotting
  p <- ggplot(data_p, aes(y = TEVp_AA, x = values, fill = TEVp_AA)) +
    geom_col(aes(fill = ifelse(AAwt == "w", "wt", ifelse(TEVp_AA == variantAA, "reported \nvariant", "other")))) +
    # geom_point() +
    scale_fill_manual(values = c("reported \nvariant" = "#ca0020", 
                                 # "max" = "#ca0020",
                                 "wt" = "#404040", 
                                 "other" = "#bababa")) +
    theme_cowplot(font_size = 12) +
    facet_wrap(~TEVp_pos, nrow = ifelse(length(unique(data_p$TEVp_pos)) >= 8, 2, 1)) +
    scale_y_discrete(name = "TEVp AA") +
    scale_x_continuous(name = activity, breaks = c(0, middle_value, max_value),
                       # specify that 0 should be labeled as 0 without decimals
                       labels = function(x) ifelse(x == 0, "0", as.character(x))) +
    theme(axis.line = element_blank()) +
    labs(fill = "TEVp AA") +
    # Annotate within each facet based on TEVp_AA and dnAUC values
    geom_text(data = variant_info, aes(label = "*", 
                                       x = ifelse(is.na(dnAUC), 0, dnAUC)),
              size = 7, hjust = -0.1, vjust = 0.8)
  
  # Add a plot title
  p <- p + ggtitle((paste0(variant, "\nactivity using ", TEVsAA, " in ", TEVspos))) +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
  # p <- p + ggtitle(paste(paste0(variant, " activity using ", TEVsAA, " in ", TEVspos, "\n"), paste(variant_info$mutation, collapse = ", "))) +
  #   theme(plot.title = element_text(hjust = 0.5))
  
  # show plot
  print(p)
  # save plot to file
  save_plot(p, filename = paste0("./05_DMS/plots/literature_variants/check_", variant, "_", TEVspos, "_", TEVsAA, "_", activity, ".pdf"),
            base_height = height, 
            base_asp = asp_ratio)
}


# write a function to plot a positional specificity profile for a set of variants
spec_prof_positional <- function(
  TEVspos, TEVp_var, activity, legend, plottitle, file_name, height, asp_ratio){
  # preparing data for the plot
  plot_data <- data_heatmaps %>%
    # add TEVp controls
    full_join(., data_controls) %>%
    filter(TEVs_pos == TEVspos,
           TEVp %in% TEVp_var,
           TEVs_AA != "*") %>%
    # remove copy entries in case TEVp wt is among the variants to be plotted (i.e. in TEVp_var)
    select(-TEVp_pos, -TEVp_AA) %>%
    distinct()
    
  # add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
  # this is important to have a grey box in the heatmap for missing values
  # create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations - 3381 lines
  df_combinations <- expand.grid(
    TEVp = unique(plot_data$TEVp),
    TEVs_AA = unique(aa_properties$AA))
  # merge with data of actually measured TEVp x TEVs combinations
  plot_data <- plot_data %>%
    full_join(., df_combinations,
              by = c("TEVp", "TEVs_AA")) %>%
    # give the name "values" to the column that should be plotted and is defined by "activity"
    rename(values = all_of(activity))
  
  # plotting
  p <- ggplot(plot_data, aes(y = factor(TEVp, levels = rev(TEVp_var)), 
                             x = factor(TEVs_AA, levels = sort(unique(TEVs_AA))),
                             fill = values)) +
    geom_tile() + # make heat map
    scale_fill_viridis(option="magma", limits = c(0, max(plot_data$values))) +
    theme_cowplot(font_size = 12) +
    ylab("") +
    theme(axis.line = element_blank()) +
    labs(fill = legend) +
    scale_x_discrete(name = paste("Amino acid in", TEVspos), expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    ggtitle(plottitle) +
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    theme(axis.title.y = element_blank()) +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          # legend.key.height = unit(10, "pt"),
          legend.key.width = unit(10, "pt")) #+
    # theme(legend.position = "top", legend.justification = "center")
  
  # show plot
  print(p)
  
  # save plot to file
  save_plot(p, filename = paste0("./05_DMS/plots/literature_variants/", file_name, ".pdf"),
            base_height = height,
            base_asp = asp_ratio)
}


## plotting of literature variants ----


# plot N176T (Carrico et al. 2016)
check_plot(variant = "N176T",
           TEVspos = "P6",
           TEVsAA = "T",
           activity = "dnAUC",
           height = 3.71,
           asp_ratio = 0.8)

# plot NX1 (Pestalozzi)
for (i in c("W", "Q", "M")) {
  check_plot(variant = "NX1",
             TEVspos = "P1'",
             TEVsAA = i,
             activity = "dnAUC",
             height = 3.71,
             asp_ratio = 1.618)
}

# plot I (Pestalozzi)
for (i in c("D", "E")) {
  check_plot(variant = "I",
             TEVspos = "P6",
             TEVsAA = i,
             activity = "dnAUC",
             height = 3.71,
             asp_ratio = 1)
}


# plot positional specificity profile of N176I and S120R mutations (Pestalozzi)
spec_prof_positional(TEVp_var = c("wt", "I", 
                                  literature_variants_long %>%
                                    filter(Variant_name == "I") %>%
                                    pull(mutation)),
                     TEVspos = "P6",
                     activity = "dnAUCrel",
                     legend = "Relative \nactivity", 
                     plottitle = "Effect of TEVp I mutations",
                     file_name = "check_TEVpI_P6_rel", 
                     height = 2.5, 
                     asp_ratio = 1.6)


# plot N176T (Carrico et al. 2016)
check_plot(variant = "N176T",
           TEVspos = "P6",
           TEVsAA = "T",
           activity = "dnAUC",
           height = 3.71,
           asp_ratio = 0.8)

# plot N171D (Carrico et al. 2016)
check_plot(variant = "N171D",
           TEVspos = "P6",
           TEVsAA = "P",
           activity = "dnAUC",
           height = 3.71,
           asp_ratio = 0.8)

# plot PE10 (Yi et al. 2013)
check_plot(variant = "PE10",
           TEVspos = "P1",
           TEVsAA = "E",
           activity = "dnAUC",
           height = 3.71,
           asp_ratio = 2)

# plot positional specificity profile of PE10 (Yi et al. 2013)

spec_prof_positional(TEVp_var = c("wt",
                                  literature_variants_long %>%
                                  filter(Variant_name == "PE10") %>%
                                  pull(mutation),
                                  "V219P"),
                     TEVspos = "P1",
                     activity = "dnAUCrel",
                     legend = "Relative \nactivity",
                     plottitle = paste("Effect of PE10 mutations"),
                     file_name = paste0("check_PE10_P1_rel"),
                     height = 3.4,
                     asp_ratio = 1.2)

# plot positional P1 specificity profile of PE10 (Yi et al. 2013)
spec_prof_positional(TEVp_var = c("wt", 
                                  literature_variants_long %>%
                                    filter(Variant_name == "PE10") %>%
                                    pull(mutation),
                                  "V219P"),
                     TEVspos = "P1",
                     activity = "dnAUCrel",
                     legend = "Relative \nactivity", 
                     plottitle = "Effect of PE10 mutations",
                     file_name = "check_PE10_P1_rel", 
                     height = 3.4, 
                     asp_ratio = 1.2)


# # plot PH7 (Yi et al. 2013)
# check_plot(variant = "PH7",
#            TEVspos = "P1",
#            TEVsAA = "H",
#            activity = "dnAUC",
#            height = 3.71,
#            asp_ratio = 1.35)
# 
# # plot PE21 (Yi et al. 2013)
# check_plot(variant = "PH21",
#            TEVspos = "P1",
#            TEVsAA = "H",
#            activity = "dnAUC",
#            height = 6.8,
#            asp_ratio = 1)
# 
# # plot M3 (Verhoeven et al. 2012)
# check_plot(variant = "M3",
#            TEVspos = "P1'",
#            TEVsAA = "D",
#            activity = "dnAUC",
#            height = 6.8,
#            asp_ratio = 1)
# 
# # plot A01 (Meister et al. 2023)
# check_plot(variant = "A01",
#            TEVspos = "P3",
#            TEVsAA = "V",
#            activity = "dnAUC",
#            height = 3.71,
#            asp_ratio = 2)
# 
# # plot C04 (Meister et al. 2023)
# check_plot(variant = "C03",
#            TEVspos = "P3",
#            TEVsAA = "V",
#            activity = "dnAUC",
#            height = 3.71,
#            asp_ratio = 2)
# 
# # plot E05 (Meister et al. 2023)
# check_plot(variant = "E05",
#            TEVspos = "P3",
#            TEVsAA = "V",
#            activity = "dnAUC",
#            height = 3.71,
#            asp_ratio = 2)
# 
# # plot E05 (Meister et al. 2023)
# check_plot(variant = "E05",
#            TEVspos = "P5",
#            TEVsAA = "K",
#            activity = "dnAUC",
#            height = 3.71,
#            asp_ratio = 2)
# 
# # plot E05 (Meister et al. 2023)
# check_plot(variant = "E05",
#            TEVspos = "P1'",
#            TEVsAA = "A",
#            activity = "dnAUC",
#            height = 3.71,
#            asp_ratio = 2)

# check variant P6 (Packer et al. 2017)
spec_prof_positional(TEVp_var = c("wt", 
                                  literature_variants_long %>%
                                    filter(Variant_name == "P6") %>%
                                    pull(mutation)),
                     TEVspos = "P6",
                     activity = "dnAUCrel",
                     legend = "Relative \nactivity", 
                     plottitle = "Individual mutations of variant P6",
                     file_name = paste0("check_P6_P6_rel"), 
                     height = 2.5, 
                     asp_ratio = 1.6)


# plot positional specificity profile of all TEVp NX1 mutations (Pestalozzi)
spec_prof_positional(TEVp_var = c("VI", "NX1", 
                                  literature_variants_long %>%
                                    filter(Variant_name == "NX1") %>%
                                    pull(mutation)),
                     TEVspos = "P1'",
                     activity = "dnAUCrel",
                     legend = "Relative \nactivity", 
                     plottitle = "Effect of individual mutations of variant NX1",
                     file_name = "check_NX1_P1'_rel", 
                     height = 3.55, 
                     asp_ratio = 1.301408451)

# plot positional specificity profile of all TEVp NX1 mutations (Pestalozzi)
spec_prof_positional(TEVp_var = c("VI", "NX1", 
                                  literature_variants_long %>%
                                    filter(Variant_name == "NX1") %>%
                                    pull(mutation)),
                     TEVspos = "P3",
                     activity = "dnAUCrel",
                     legend = "Relative \nactivity", 
                     plottitle = "Effect of individual mutations of variant NX1",
                     file_name = "check_NX1_P3_rel", 
                     height = 2.1, 
                     asp_ratio = 2.2)


# plot P5 specificity profile of truncated mutants
spec_prof_positional(TEVp_var = c("wt", "V219*" 
                                  # data %>%
                                  #   arrange(TEVp_pos) %>%
                                  #   filter(!TEVp_pos %in% c(220, 222, 232, 233)) %>%
                                  #   filter(TEVp_pos >= 219,
                                  #          TEVp_AA == "*") %>%
                                  #   pull(TEVp) %>% unique()
                                  ),
                     TEVspos = "P5",
                     activity = "dnAUCrel",
                     legend = "Relative \nactivity", 
                     plottitle = "Effect of truncating TEVp at the C-terminus",
                     file_name = "truncated_variants_P5_rel", 
                     height = 2.1, 
                     asp_ratio = 2.2)

# Getting the top hits ----------------------------------------------------

### Top TEVp positions (avg) ----

# define an exemplary  substrate of interest
TEVspos <- "P1'"
TEVsAA <- "W"

# define a cutoff for the AUC
cutoff <- 0.15

# define a threshold for the number of mutants tested per position
threshold <- 5

# summarize the data - this is just to get an overview
summary <- data_heatmaps %>%
  filter(TEVs_pos == TEVspos,
         TEVs_AA == TEVsAA,
         AUC <= cutoff,
         TEVp_AA != "*",
         TEVp_pos %in% seq(1, 234)) %>%
  group_by(TEVp_pos) %>%
  summarize(avg_dnAUC = mean(dnAUC),
            min = min(dnAUC),
            max = max(dnAUC),
            AAmax = TEVp_AA[which.max(dnAUC)],
            obs = n(),
            obs_rel = n()/20)

# get the data of interest
exempl <- data_heatmaps %>%
  filter(TEVs_pos == TEVspos,
         TEVs_AA == TEVsAA,
         TEVp_AA != "*") %>%
  # add information which TEVp_AA yields the highest dnAUC
  group_by(TEVp_pos) %>%
  mutate(AAmax = TEVp_AA[which.max(dnAUC)]) %>%
  ungroup

# define the top TEVp positions (highest average activity)
TEVppos <- exempl %>%
  group_by(TEVp_pos) %>%
  summarize(avg_dnAUC = mean(dnAUC),
            obs = n()) %>%
  filter(obs >= threshold) %>%
  arrange(desc(avg_dnAUC)) %>%
  head(5) %>% 
  pull(TEVp_pos)

# filter for top positions
exempl <- filter(exempl, TEVp_pos %in% TEVppos)

# give the inactive cases a very small activity, so they can be 
# distinguished from cases that were not tested
exempl$dnAUC[exempl$dnAUC == 0] <- 0.002

# add the wt AA to the data frame to color it differently later
exempl <- left_join(exempl, 
                    TEVpwt %>% 
                      select(-TEVp) %>% 
                      rename(AAwt = TEVp_AA), 
                    by = "TEVp_pos")

# plotting
p <- ggplot(exempl, aes(y = TEVp_AA, x = dnAUC, fill = TEVp_AA)) +
  geom_col(aes(fill = ifelse(TEVp_AA == AAmax, "max", ifelse(TEVp_AA == AAwt, "wt", "other")))) +
  scale_fill_manual(values = c("max" = "#ca0020", "wt" = "#404040", "other" = "#bababa")) +
  theme_cowplot(font_size = 12) +
  facet_wrap(~factor(TEVp_pos, levels = TEVppos), ncol = 5) +
  scale_y_discrete(name = "TEVp AA") +
  theme(axis.line = element_blank()) +
  labs(fill = "TEVp AA") +
  theme(axis.text.x = element_text(size = 8)) + # change size of x axis labels
  ggtitle(paste0("Top TEVp positions for ", TEVsAA, " in ", TEVspos)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(data = exempl %>% group_by(TEVp_pos) %>% summarize(avg_dnAUC = mean(dnAUC)),
             aes(xintercept = avg_dnAUC),
             linetype = "dashed", color = "black", size = 0.5)

# show plot
print(p)


#### Corresponding spec. prof. ----

# define the top hits

# define the top TEVp positions (highest average activity)
hits <- exempl %>%
  group_by(TEVp_pos) %>%
  summarize(avg_dnAUC = mean(dnAUC),
            TEVp = TEVp[which.max(dnAUC)]) %>%
  arrange(desc(avg_dnAUC)) %>%
  pull(TEVp)

# get the data of the top hits
specificity <- data_heatmaps %>%
  filter(TEVp %in% hits,
         TEVs_AA != "*")

# add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
# this is important to have a grey box in the heatmap for missing values
# create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations
df_combinations <- expand.grid(
  # TEVp_AA = unique(specificity$TEVp_AA),
  TEVs_AA = unique(specificity$TEVs_AA),
  # TEVp_pos = unique(specificity$TEVp_pos),
  TEVp = hits,
  TEVs_pos = unique(specificity$TEVs_pos))
# merge with data of actually measured TEVp x TEVs combinations
specificity <- specificity %>%
  full_join(., df_combinations)

# specificity <- filter(specificity, !is.na(TEVp_pos))

# plotting
p2 <- ggplot(specificity, aes(x = factor(TEVs_pos, levels = rev(c("P1'", "P1", "P2", "P3", "P4", "P5", "P6"))),
                              y = TEVs_AA, fill = dnAUC)) +
  geom_tile() +
  facet_wrap(~TEVp_pos, ncol = 5) +
  facet_wrap(~factor(TEVp, levels = hits), ncol = 5) +
  scale_fill_viridis(option="magma") +
  theme_cowplot(font_size = 12) +
  scale_y_discrete(name = "TEVs AA", expand = c(0, 0)) +
  scale_x_discrete(name = "TEVs position", expand = c(0, 0))+
  theme(axis.line = element_blank()) +
  theme(axis.text.x = element_text(size = 8)) # change size of x axis labels

# show plot
print(p2)

# put both plots underneath each other
plotgrid <- plot_grid(p, p2, ncol = 1, labels = "auto")
# show grid
print(plotgrid)

# save plot to file
save_plot(plotgrid, filename = paste0("./05_DMS/plots/top_positions_", TEVspos, TEVsAA, "_grid.pdf"),
          base_height = 7.42, 
          base_asp = 1.2)



### Top TEVp mutants ----

# define an exemplary  substrate of interest
TEVspos <- "P1'"
TEVsAA <- "W"

# define a cutoff for the AUC
cutoff <- 0.15

# get the data of interest
exempl <- data_heatmaps %>%
  filter(TEVs_pos == TEVspos,
         TEVs_AA == TEVsAA,
         AUC <= cutoff,
         TEVp_AA != "*") %>%
  # add information which TEVp_AA yields the highest dnAUC
  group_by(TEVp_pos) %>%
  mutate(AAmax = TEVp_AA[which.max(dnAUC)]) %>%
  ungroup

# define the TEVp positions with top TEVps (highest single variant)
TEVppos <- exempl %>%
  arrange(desc(dnAUC)) %>%
  pull(TEVp_pos) %>%
  unique() %>%
  .[seq(1, 5)] # take the 5 highest

# filter for top positions
exempl <- filter(exempl, TEVp_pos %in% TEVppos)

# give the inactive cases a very small activity, so they can be
# distinguished from cases that were not tested
exempl$dnAUC[exempl$dnAUC == 0] <- 0.002

# add the wt AA to the data frame to color it differently later in the plot
exempl <- left_join(exempl, 
                    TEVpwt %>% 
                      select(-TEVp) %>% 
                      rename(AAwt = TEVp_AA), 
                    by = "TEVp_pos")

# plotting
p <- ggplot(exempl, aes(y = TEVp_AA, x = dnAUC, fill = TEVp_AA)) +
  geom_col(aes(fill = ifelse(TEVp_AA == AAmax, "max", ifelse(TEVp_AA == AAwt, "wt", "other")))) +
  scale_fill_manual(values = c("max" = "#ca0020", "wt" = "#404040", "other" = "#bababa")) +
  theme_cowplot(font_size = 12) +
  facet_wrap(~factor(TEVp_pos, levels = TEVppos), ncol = 5) +
  scale_y_discrete(name = "TEVp AA") +
  theme(axis.line = element_blank()) +
  labs(fill = "TEVp AA") +
  theme(axis.text.x = element_text(size = 8)) + # change size of x axis labels
  ggtitle(paste0("Top TEVp mutants for ", TEVsAA, " in ", TEVspos)) +
  theme(plot.title = element_text(hjust = 0.5))


# show plot
print(p)
# save plot to file
save_plot(p, filename = paste0("./05_DMS/plots/top_mutants_", TEVspos, TEVsAA, ".pdf"),
          base_height = 3.71, 
          base_asp = 1.618)

#### Corresponding spec. prof. ----

# define the top hits
hits <- exempl %>%
  filter(TEVp_pos %in% TEVppos,
         TEVp_AA == AAmax) %>%
  arrange(desc(dnAUC)) %>%
  pull(TEVp)

# get the data of the top hits
specificity <- data_heatmaps %>%
  filter(TEVp %in% hits,
         TEVs_AA != "*")

# add the TEVp x TEVs combinations that were not covered in the screen and set the value to NA
# this is important to have a grey box in the heatmap for missing values
# create mock data frame with all unique TEVp, TEVp position and TEVs AA combinations
df_combinations <- expand.grid(
  # TEVp_AA = unique(specificity$TEVp_AA),
  TEVs_AA = unique(specificity$TEVs_AA),
  # TEVp_pos = unique(specificity$TEVp_pos),
  TEVp = hits,
  TEVs_pos = unique(specificity$TEVs_pos))
# merge with data of actually measured TEVp x TEVs combinations
specificity <- specificity %>%
  full_join(., df_combinations)

# specificity <- filter(specificity, !is.na(TEVp_pos))

# plotting
p2 <- ggplot(specificity, aes(x = factor(TEVs_pos, levels = rev(c("P1'", "P1", "P2", "P3", "P4", "P5", "P6"))),
                             y = TEVs_AA, fill = dnAUC)) +
  geom_tile() +
  facet_wrap(~TEVp_pos, ncol = 5) +
  facet_wrap(~factor(TEVp, levels = hits), ncol = 5) +
  scale_fill_viridis(option = "magma", 
                     # limits = c(0, 0.13)
                     ) +
  theme_cowplot(font_size = 12) +
  scale_y_discrete(name = "TEVs AA", expand = c(0, 0)) +
  scale_x_discrete(name = "TEVs position", expand = c(0, 0))+
  theme(axis.line = element_blank())
# show plot
print(p2)


# put both plots underneath each other
plotgrid <- plot_grid(p, p2, ncol = 1, labels = "auto")
# show grid
print(plotgrid)

# save plot to file
save_plot(plotgrid, filename = paste0("./05_DMS/plots/top_mutants_", TEVspos, TEVsAA, "_grid.pdf"),
          base_height = 11, 
          base_asp = 0.9)








