# load libraries
library(tidyverse)
library(viridis)
library(cowplot)
library(stringi)
library(matrixStats)


# create folders ----------------------------------------------------------

# create a folder to store the individual plots and data in
dir.create("./02_DMS_controls_false_positives")
dir.create("./02_DMS_controls_false_positives/results")
dir.create("./02_DMS_controls_false_positives/plots")


# prepare data ------------------------------------------------------------

# read data
# data <- read.table("./data/Data_DMS_controls.txt", sep = "\t", header = TRUE)

data_BCs <- read.table("./data/Data_DMS_controls_unmapped.txt", sep = "\t", header = TRUE)

# get mutated TEVs BCs (i.e., mutated backbones) identified in another analysis
mutated_backbones <- read.table("./01.1_DMS_mutated_backbones/results/mutated_backbones.txt", sep = "\t", header = TRUE)

# throw out mutated backbones 
data_BCs <- filter(data_BCs, !BC_TEVs %in% mutated_backbones$BC_TEVs)

# throw out columns not needed
data_BCs <- data_BCs %>%
  select(-c(flask, source_TEVp, list_TEVs, reason_TEVs, source_TEVs, libraryID, PlasmidID))

# make a min column specifying the minimum read per time point
data_BCs$min <- rowMins(as.matrix(data_BCs[, c(paste0('N', 1:5))] + data_BCs[, c(paste0('F', 1:5))]))

### calculate fraction flipped and AUC  ----

# calculate fraction flipped for each TEVp BC - TEVs BC combination
df_fraction <- data_BCs[, c(paste0('F', 1:5))] / (data_BCs[, c(paste0('N', 1:5))] + data_BCs[, c(paste0('F', 1:5))])
names(df_fraction) <- paste0('t', 1:5)

# add the fraction flipped columns to the data
data_BCs$t1 <- round(df_fraction$t1, 3)
data_BCs$t2 <- round(df_fraction$t2, 3)
data_BCs$t3 <- round(df_fraction$t3, 3)
data_BCs$t4 <- round(df_fraction$t4, 3)
data_BCs$t5 <- round(df_fraction$t5, 3)

### calculate AUC
x <- as.integer(0:4) # get time points

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_fraction), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_fraction[,i]) + (x[i+1]-x[i])*(df_fraction[,i+1]-df_fraction[,i])*0.5
  m_auc[,i] <- auc
}

# calculate AUCs and add columns to data
data_BCs$AUC <- round(rowSums(m_auc[,1:4])/4, 3)



# Numbers before removing outliers ---------------------------------------------

# count how many unique TEVp-TEVs combinations we covered with min 20 reads on protein level
number_min20protein <- data %>%
  filter(n_mut_TEVs <= 1) %>%
  mutate(TEVp_TEVs = paste0(TEVp, "_", TEVs)) %>%
  pull(TEVp_TEVs) %>%
  unique() %>%
  length()

cat("Number of unique TEVp-TEVs combinations (single mutants only) applying a min read threshold of 20 per time point on protein level:",
    number_min20protein)


# count how many unique TEVp-TEVs combinations we covered with min 20 reads on DNA level
number_min20DNA <- data_BCs %>%
  filter(n_mut_TEVs <= 1,
         min >= 20) %>%
  mutate(TEVp_TEVs = paste0(TEVp, "_", TEVs)) %>%
  pull(TEVp_TEVs) %>%
  unique() %>%
  length()

cat("Number of unique TEVp-TEVs combinations (single mutants only) applying a min read threshold of 20 per time point on DNA level:",
    number_min20DNA)

ratio <- 100 - number_min20DNA / number_min20protein *100
# Print the result
cat("Percent of unique TEVp-TEVs combinations lost when applying a min 20 read threshold per time point for each plasmid:",
    ratio)

# Check for potentially mutated TEVs backbones -------------------------------

# Search for TEVs BCs (i.e. individual plasmid backbones) that consistently yield
# high AUCs irrespective of the TEVp variant

## inspection of TEVs BCs ----

# filter for min 20 reads and exclude the wt TEVs BC, i.e. the plasmid with wt TEVs
# that was used to clone the TEVp library
my_data <- data_BCs%>%
  filter(min >= 20) %>%
  select(BC_TEVp, BC_TEVs, AUC) #%>%
  # throw out wt TEVs BC
  # filter(BC_TEVs != "ACGTCGCTGA")

# calculate the median AUC for each TEVs_BC
TEVs_BC_medians <- my_data %>%
  group_by(BC_TEVs) %>%
  summarize(median = median(AUC),
            n = n()) %>%
  arrange(desc(median)) %>%
  # add TEVs column
  left_join(., select(data_BCs, TEVs, BC_TEVs) %>% distinct(), by = "BC_TEVs") %>%
  # throw out all substrates containing a stop codon
  filter(!grepl("\\*", TEVs)) %>%
  mutate(rank = row_number())


# add the count (n) and TEVs column to the data for plotting reasons later
my_data <- left_join(my_data,
                     TEVs_BC_medians %>% select(BC_TEVs, n), by = "BC_TEVs") %>%
  # add TEVs column
  left_join(., select(data_BCs, TEVs, BC_TEVs) %>% distinct(), by = "BC_TEVs")

### highest TEVs_BCs ----

# manually define how many TEVs BCs to plot and order them according to their median
TEVs_BC_order <- TEVs_BC_medians %>%
  pull(BC_TEVs) %>%
  .[1:50]

# plot the TEVs plasmid backbones
p <- ggplot(data = my_data,
            # data = my_data %>% filter(BC_TEVs %in% TEVs_BC_medians$BC_TEVs[1:50]),
            aes(x = factor(BC_TEVs, levels=TEVs_BC_order), y = AUC, fill = n)) +
  # geom_violin() +
  geom_boxplot(width = 0.8, color = "black", alpha = 1) +
  labs(x = "TEVs BC (= plasmid backbone ID)", y = "AUC") +
  theme_cowplot(font_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# print plot
print(p)


# They all look good, none are "off".


### lowest TEVs_BCs ----

# do the same with very low TEVs_BCs...

# get lowest TEVs BCs
TEVs_BC_order <- TEVs_BC_medians %>%
  arrange(median) %>%
  pull(BC_TEVs) %>%
  .[1:50]

# plot the TEVs plasmid backbones
p <- ggplot(data = my_data,
            # data = my_data %>% filter(BC_TEVs %in% TEVs_BC_medians$BC_TEVs[1:50]),
            aes(x = factor(BC_TEVs, levels=TEVs_BC_order), y = AUC, fill = n)) +
  # geom_violin() +
  geom_boxplot(width = 0.8, color = "black", alpha = 1) +
  labs(x = "TEVs BC (= plasmid backbone ID)", y = "AUC") +
  theme_cowplot(font_size = 10) +
  coord_cartesian(ylim = c(0, 0.25)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# print plot
print(p)


# They all look good, no TEVs BC is really 'off'.



# Calculate AUCs after merging BCs per TEVp variant -----------------------


# calculate the sum of N and F per time point for each TEVp-TEVs combination
data <- aggregate(cbind(N1, N2, N3, N4, N5, F1, F2, F3, F4, F5) ~ TEVp + TEVs, data = data_BCs, FUN = sum)

# make fraction
df_fraction <- data[, c(paste0('F', 1:5))] / (data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))])
names(df_fraction) <- paste0('t', 1:5)

# add the fraction flipped columns to the data
data$t1 <- round(df_fraction$t1, 3)
data$t2 <- round(df_fraction$t2, 3)
data$t3 <- round(df_fraction$t3, 3)
data$t4 <- round(df_fraction$t4, 3)
data$t5 <- round(df_fraction$t5, 3)

### calculate AUC
x <- as.integer(0:4) # get time points

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_fraction), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_fraction[,i]) + (x[i+1]-x[i])*(df_fraction[,i+1]-df_fraction[,i])*0.5
  m_auc[,i] <- auc
}

# calculate AUC and add column to data
data$AUC <- round(rowSums(m_auc[,1:4])/4, 4)

# make a min column specifying the minimum read per time point
data$min <- rowMins(as.matrix(data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))]))

# make a sum column specifying the sum of reads over all time points
data$sum <- rowSums(data[, (c(paste0('N', 1:5), paste0('F', 1:5)))])

# throw out everything with less than 20 reads per time point
data <- filter(data, min >= 20) %>%
  arrange(desc(AUC))

# add important columns
data <- inner_join(data,
                        data_BCs %>% select(TEVp, TEVs, n_mut_TEVs) %>% distinct())


# Histogram before outlier correction ----

## min20protein ----

# count the number of replicates per TEVp-TEVs combination
replicate_samples <- data_BCs %>%
  # for this plot, samples should be considered replicates only if their coverage
  # was high enough to calculate an AUC
  filter(is.finite(AUC)) %>%
  mutate(TEVp_TEVs = paste0(TEVp, "_", TEVs)) %>%
  group_by(TEVp, TEVs) %>% 
  summarize(n_replicates = n()) %>%
  # filter out wt TEVp and wt TEVs because they are extremely overrepresented
  filter(TEVp != "wt",
         TEVs != "ENLYFQS") %>%
  arrange(desc(n_replicates))

# make a histogram plot of replicates
p <- ggplot(replicate_samples, aes(x = n_replicates)) +
  geom_histogram(binwidth = 1, fill = "#8da0cb", color = "black", alpha = 0.7) +
  scale_x_continuous(breaks = unique(replicate_samples$n_replicates)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_cowplot(font_size = 12) +
  labs(x = "Number of replicates", y = "Count of unique TEVp-TEVs combinations") +
  geom_text(stat = 'count', aes(label = ..count..), size = 3, angle = 90, vjust = 0.5, hjust = -0.2)
# print plot  
print(p)

#save plot
save_plot(p, filename = "./02_DMS_controls_false_positives/plots/histogram_replicates_before_removing_outliers_min20protein.pdf",
          base_height = 3.71, 
          base_asp = 1)

## min20DNA ----

# count the number of replicates per TEVp-TEVs combination
replicate_samples_min20DNA <- data_BCs %>%
  # for this plot, samples should be considered replicates only if their coverage
  # was high enough to calculate an AUC
  filter(is.finite(AUC)) %>%
  filter(min >= 20) %>%
  mutate(TEVp_TEVs = paste0(TEVp, "_", TEVs)) %>%
  group_by(TEVp, TEVs) %>% 
  summarize(n_replicates = n()) %>%
  # filter out wt TEVp and wt TEVs because they are extremely overrepresented
  filter(TEVp != "wt",
         TEVs != "ENLYFQS") %>%
  arrange(desc(n_replicates))

# make a histogram plot of replicates when min20 threshold is applied
p <- ggplot(replicate_samples_min20DNA, aes(x = n_replicates)) +
  geom_histogram(binwidth = 1, fill = "#8da0cb", color = "black", alpha = 0.7) +
  scale_x_continuous(breaks = unique(replicate_samples_min20DNA$n_replicates)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_cowplot(font_size = 12) +
  labs(x = "Number of replicates", y = "Count of unique TEVp-TEVs combinations") +
  geom_text(stat = 'count', aes(label = ..count..), size = 3, angle = 90, vjust = 0.5, hjust = -0.2)
# print plot  
print(p)

#save plot
save_plot(p, filename = "./02_DMS_controls_false_positives/plots/histogram_replicates_before_removing_outliers_min20DNA.pdf",
          base_height = 3.71, 
          base_asp = 1)



### Now, let's check for potential false positives, i.e. individual flipping 
### curves, that yield much higher AUCs compared to the remaining curves of the 
### same TEVp-TEVs combination.

# Data prep Option 1: min20protein ---------------------------------------------

# require minimum 20 reads on protein level (i.e. per TEVp - TEVs combination)

# Using all replicates, even when below 20 reads per time point, and 
# applying a threshold fc of e.g. 2 meaning that the top flipping curve must not
# increase the read adjusted average by more than 2-fold.
# One criteria is of course, that after removing the top curves, the remaining 
# curves must still reach a minimum of 20 reads per time point when taken together.

## deal with duplicates --------------------------------------

# whenever a sample is measured in duplicates and the two duplicates are highly
# different from each other, we cannot tell which replicate of the two is more
# correct. In that case, the sample should be removed overall.

# extract only TEVp-TEVs combinations that are represented by 2 BCs
duplicates_min20protein <- data_BCs %>%
  # only consider samples where the coverage was high enough to calculate an AUC
  filter(is.finite(AUC)) %>%
  # within TEVp-TEVs combinations
  group_by(TEVp, TEVs) %>% 
  # at exact 2 entries must exist, i.e. 2 different plasmids 
  filter(n() == 2)

# calculate the sum of N and F per time point for each TEVp-TEVs combination
df_sumNF <- aggregate(cbind(N1, N2, N3, N4, N5, F1, F2, F3, F4, F5) ~ TEVp + TEVs, data = duplicates_min20protein, FUN = sum) 

# Filter by counts: 20 reads on protein level
df_sumNF <- df_sumNF %>% filter(N1+F1 >= 20 & N2+F2 >= 20 & N3+F3 >= 20 & N4+F4 >= 20 & N5+F5 >= 20)

# Merge with original data frame to add BC_TEVp column back
duplicates_min20protein <- inner_join(df_sumNF %>% select(TEVp, TEVs),
                                      duplicates_min20protein %>% select(TEVp, TEVs, BC_TEVp, BC_TEVs, paste0("t", 1:5), paste0("N", 1:5), paste0("F", 1:5), sum, min, AUC),
                                      by = c("TEVp", "TEVs")) %>% # add a column specifying the sample
  mutate(TEVp_TEVs = paste0(TEVp, "_", TEVs))

# Next, very diverging duplicates should be detected. For this, we calculate how 
# much the upper curve differs from the lower curve. However, we only want do this
# with curves that we can calculate a reasonable AUC from, i.e. curves with at least
# 20 reads per time point. 
# For duplicates below 20 reads per time point, we decided to keep them and
# calculate the AUC over both duplicates. 

# calculate the fc of the upper curve to the lower curve
duplicates_min20protein_20 <- duplicates_min20protein %>%
  # filter only for high coverage variants! Reason is described above
  filter(min >= 20) %>%
  group_by(TEVp, TEVs) %>%
  # should still be present in duplicates
  filter(n() == 2) %>%
  mutate(fc = max(AUC)/min(AUC))

# Now, we want to identify the highly diverging duplicates based on the AUC fold
# change of the upper curve compared to the lower curve. Since the fold change
# tends to be very high when dividing small numbers, we only want to consider
# samples where at least one replicate is above the noise cutoff.

# get information about noise (this is the noise of P1' controls only)
noise <- read.table("./01_internal_controls/results/df_noise.txt", sep = "\t", header = TRUE)
cutoff <- noise %>% pull(cutoff)

# define a threshold for divergence, i.e. fold change of upper curve to lower curve
threshold <- 3

# identify samples with diverging duplicates
diverging_duplicates_min20protein <- duplicates_min20protein_20 %>%
  # consider only cases where at least one duplicate is above the noise
  filter(AUC >= cutoff, 
         fc >= threshold) %>%
  pull(TEVp_TEVs) %>%
  unique()

cat("Number of duplicates on protein level:",
    duplicates_min20protein$TEVp_TEVs %>% unique %>% length())

cat("Number of duplicates on protein level with >20 reads per replicate:",
    duplicates_min20protein_20$TEVp_TEVs %>% unique %>% length())

cat("Number of duplicates that need to be removed because they diverge a lot from each other:",
    diverging_duplicates_min20protein %>% length())

cat("Percent of duplicates that need to be removed because they diverge a lot from each other:",
    diverging_duplicates_min20protein %>% length() / duplicates_min20protein_20$TEVp_TEVs %>% unique %>% length() * 100)

## deal with replicates (n >= 3) ----

# do this analysis only on active samples, i.e. TEVp-TEVs combination above the noise cutoff

# get information about noise (this is the noise of P1' controls only)
noise <- read.table("./01_internal_controls/results/df_noise.txt", sep = "\t", header = TRUE)
cutoff <- noise %>% pull(cutoff)

# get a table of all active variants
active <- filter(data, AUC >= cutoff) %>%
  select(TEVp, TEVs)

# sort out all inactive TEVp-TEVs combinations from the data_BCs table
data_BCs_active <- inner_join(active, data_BCs, by = c("TEVp", "TEVs"))

# extract only TEVp-TEVs combinations that are represented with at least 3 BCs
over3BCs <- data_BCs_active %>%
  # only consider samples where the coverage was high enough to calculate an AUC
  filter(is.finite(AUC)) %>%
  # within TEVp-TEVs combinations
  group_by(TEVp, TEVs) %>% 
  # at least 3 entries must exist, i.e. 3 different plasmids
  filter(n() >= 3) 

## Filter for TEVp-TEVs combinations that have at least 20 reads per time point (on protein level!)

# calculate the sum of N and F per time point for each TEVp-TEVs combination
df_sumNF <- aggregate(cbind(N1, N2, N3, N4, N5, F1, F2, F3, F4, F5) ~ TEVp + TEVs, data = over3BCs, FUN = sum) 

# Filter by counts: 20 reads on protein level
df_sumNF <- df_sumNF %>% filter(N1+F1 >= 20 & N2+F2 >= 20 & N3+F3 >= 20 & N4+F4 >= 20 & N5+F5 >= 20)

# Merge with original data frame to add BC_TEVp column back
over3BCs_min20protein <- inner_join(df_sumNF %>% select(TEVp, TEVs),
                                    over3BCs %>% select(TEVp, TEVs, BC_TEVp, BC_TEVs, paste0("t", 1:5), paste0("N", 1:5), paste0("F", 1:5), sum, min, AUC), 
                                    by = c("TEVp", "TEVs")) %>%
  mutate(plasmidID = paste0(BC_TEVp, "_", BC_TEVs))



# Data prep Option 2: min20DNA -------------------------------------------------

# require minimum 20 reads on DNA level (i.e. per TEVp BC - TEVs BC combination)

# Using only replicates with minimum 20 reads per time point and 
# applying a threshold fc of e.g. 2 meaning that the top flipping curve must not
# increase the read adjusted average by more than 2-fold.

## filter for min 20 reads per time point and BC ----

data_BCs_min20DNA <- filter(data_BCs, min >= 20)

## recalculate the AUCs after filtering for min 20 reads per time point and BC ----

# calculate the sum of N and F per time point for each TEVp-TEVs combination
data_min20DNA <- aggregate(cbind(N1, N2, N3, N4, N5, F1, F2, F3, F4, F5) ~ TEVp + TEVs, data = data_BCs_min20DNA, FUN = sum) 

# make fraction
df_fraction <- data_min20DNA[, c(paste0('F', 1:5))] / (data_min20DNA[, c(paste0('N', 1:5))] + data_min20DNA[, c(paste0('F', 1:5))])
names(df_fraction) <- paste0('t', 1:5)

# add the fraction flipped columns to the data
data_min20DNA$t1 <- round(df_fraction$t1, 3)
data_min20DNA$t2 <- round(df_fraction$t2, 3)
data_min20DNA$t3 <- round(df_fraction$t3, 3)
data_min20DNA$t4 <- round(df_fraction$t4, 3)
data_min20DNA$t5 <- round(df_fraction$t5, 3)

### calculate AUC
x <- as.integer(0:4) # get time points

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_fraction), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_fraction[,i]) + (x[i+1]-x[i])*(df_fraction[,i+1]-df_fraction[,i])*0.5
  m_auc[,i] <- auc
}

# calculate AUC and add column to data
data_min20DNA$AUC <- round(rowSums(m_auc[,1:4])/4, 4)

# make a min column specifying the minimum read per time point
data_min20DNA$min <- rowMins(as.matrix(data_min20DNA[, c(paste0('N', 1:5))] + data_min20DNA[, c(paste0('F', 1:5))]))

# make a sum column specifying the sum of reads over all time points
data_min20DNA$sum <- rowSums(data_min20DNA[, (c(paste0('N', 1:5), paste0('F', 1:5)))])

# throw out everything with less than 20 reads per time point
data_min20DNA <- filter(data_min20DNA, min >= 20) %>%
  arrange(desc(AUC))

# add important columns
data_min20DNA <- inner_join(data_min20DNA,
                            data_BCs_min20DNA %>% select(TEVp, TEVs, n_mut_TEVs) %>% distinct())


## deal with duplicates --------------------------------------

# whenever a sample is measured in duplicates and the two duplicates are highly
# different from each other, we cannot tell which replicate of the two is more
# correct. In that case, the sample should be removed overall.

# extract only TEVp-TEVs combinations that are represented by 2 BCs
duplicates_min20DNA <- data_BCs_min20DNA %>%
  # only consider samples where the coverage was high enough to calculate an AUC
  filter(is.finite(AUC)) %>%
  # within TEVp-TEVs combinations
  group_by(TEVp, TEVs) %>% 
  # at exact 2 entries must exist, i.e. 2 different plasmids 
  filter(n() == 2)

# Next, very diverging duplicates should be detected. For this, we calculate how 
# much the upper curve differs from the lower curve. 

# calculate the fc of the upper curve to the lower curve
duplicates_min20DNA <- duplicates_min20DNA %>% 
  # add a column specifying the sample
  mutate(TEVp_TEVs = paste0(TEVp, "_", TEVs)) %>%
  group_by(TEVp_TEVs) %>%
  mutate(fc = max(AUC)/min(AUC))

# Now, we want to identify the highly diverging duplicates based on the AUC fold
# change of the upper curve compared to the lower curve. Since the fold change
# tends to be very high when dividing small numbers, we only want to consider
# samples where at least one replicate is above the noise cutoff.

# get information about noise (this is the noise of P1' controls only)
noise <- read.table("./01_internal_controls/results/df_noise.txt", sep = "\t", header = TRUE)
cutoff <- noise %>% pull(cutoff)

# define a threshold for divergence, i.e. fold change of upper curve to lower curve
threshold <- 3

# identify samples with diverging duplicates
diverging_duplicates_min20DNA <- duplicates_min20DNA %>%
  # consider only cases where at least one duplicate is above the noise
  filter(AUC >= cutoff, 
         fc >= threshold) %>%
  pull(TEVp_TEVs) %>%
  unique()


cat("Number of duplicates on DNA level with >20 reads per replicate:",
    duplicates_min20DNA$TEVp_TEVs %>% unique %>% length())

cat("Number of duplicates that need to be removed because they diverge a lot from each other:",
    diverging_duplicates_min20DNA %>% length())

cat("Percent of duplicates that need to be removed because they diverge a lot from each other:",
    diverging_duplicates_min20DNA %>% length() / duplicates_min20DNA$TEVp_TEVs %>% unique %>% length() * 100)



## deal with replicates (n >= 3) ----

# do this analysis only on active samples, i.e. TEVp-TEVs combination above the noise cutoff

# get information about noise (this is the noise of P1' controls only)
noise <- read.table("./01_internal_controls/results/df_noise.txt", sep = "\t", header = TRUE)
cutoff <- noise %>% pull(cutoff)

# get a table of all active variants
active <- filter(data_min20DNA, AUC >= cutoff) %>%
  select(TEVp, TEVs)

# sort out all inactive TEVp-TEVs combinations from the data_BCs table
data_BCs_active <- inner_join(active, data_BCs, by = c("TEVp", "TEVs"))

# require every biological replicate to have at least 20 reads per time point
over3BCs_min20DNA <- data_BCs_active %>%
  # every replicate on DNA level must have at least 20 reads per time point!
  filter(min >= 20) %>%
  # within TEVp-TEVs combinations
  group_by(TEVp, TEVs) %>% 
  # at least 3 entries must exist, i.e. 3 different TEVp BCs
  filter(n() >= 3) %>%
  mutate(plasmidID = paste0(BC_TEVp, "_", BC_TEVs))



# Prepare for looping through both alternatives ---------------------------

# define names for the two alternatives
alternatives <- c("min20proteinlevel", "min20DNAlevel")

# Since the underlying data differ based on the option, store the respective
# data sets in a list
list_data <- list(min20proteinlevel = data, min20DNAlevel = data_min20DNA)
list_data_BCs <- list(min20proteinlevel = data_BCs, min20DNAlevel = data_BCs_min20DNA)
list_diverging_duplicates <- list(min20proteinlevel = diverging_duplicates_min20protein, min20DNAlevel = diverging_duplicates_min20DNA)
list_biol_reps <- list(min20proteinlevel = over3BCs_min20protein, min20DNAlevel = over3BCs_min20DNA)


# Replicates (n >= 3) -----------------------------------------------------


## Calculate fc caused by top curve -------------------------------------------

# do the following steps for both alternatives, min20protein and min20DNA

# initialize a list to store the results in
list_compareAUCs <- list()

# Loop through the two alternative data sets
for (alternative in c(1:2)) {
  
  # get the correct data sets for the respective alternative
  data <- list_data[[alternative]]
  data_BCs <- list_data_BCs[[alternative]]
  biol_reps <- list_biol_reps[[alternative]]
  
  # create folder to store the plots in
  dir.create(paste0("./02_DMS_controls_false_positives/plots/", alternatives[alternative]))
  
  # Count the number of unique TEVp-TEVs combinations in this data frame
  num_unique_combinations <- biol_reps %>%
    select(TEVp, TEVs) %>%
    distinct() %>%
    nrow()
  
  # Print the result
  message <- "Number of active TEVp-TEVs combinations with at least 3 different replicates in that combination:"
  cat(message, num_unique_combinations)
  
  # save this result to a file
  write.table(data.frame(n = num_unique_combinations,
                         criteria = message),
              file = paste0("./02_DMS_controls_false_positives/results/n_biol_reps_", alternatives[alternative], ".txt"), sep = "\t", row.names=FALSE, quote = F)
  
  ### calculate AUCs after removing the top flipping curve ----
  
  # identify the plasmid with the highest AUC within each TEVp-TEVs combination
  compareAUCs <- biol_reps %>%
    group_by(TEVp, TEVs) %>%
    summarize(AUCmax = max(AUC, na.rm = TRUE),
              plasmidID_AUCmax = plasmidID[which.max(AUC)],
              replicates = n())
  
  # add the original (read adjusted) AUC
  compareAUCs <- left_join(compareAUCs,
                           data %>% select(TEVp, TEVs, AUC),
                           by = c("TEVp", "TEVs"))
  
  # remove the plasmid yielding the highest AUC per TEVp-TEVs combination
  biol_reps_nomax <- anti_join(biol_reps,
                               compareAUCs %>% rename("plasmidID" = "plasmidID_AUCmax"),
                               by = c("TEVp", "TEVs", "plasmidID"))
  
  # calculate the sum of N and F per time point for each TEVp-TEVs combination
  data_nomax <- aggregate(cbind(N1, N2, N3, N4, N5, F1, F2, F3, F4, F5) ~ TEVp + TEVs, data = biol_reps_nomax, FUN = sum)
  
  # make fraction
  df_fraction <- data_nomax[, c(paste0('F', 1:5))] / (data_nomax[, c(paste0('N', 1:5))] + data_nomax[, c(paste0('F', 1:5))])
  names(df_fraction) <- paste0('t', 1:5)
  
  # add the fraction flipped columns to the data
  data_nomax$t1 <- round(df_fraction$t1, 3)
  data_nomax$t2 <- round(df_fraction$t2, 3)
  data_nomax$t3 <- round(df_fraction$t3, 3)
  data_nomax$t4 <- round(df_fraction$t4, 3)
  data_nomax$t5 <- round(df_fraction$t5, 3)
  
  ### calculate AUC
  x <- as.integer(0:4) # get time points
  
  # initialize helper matrix for AUC calculation
  m_auc <- matrix(nrow = nrow(df_fraction), ncol = (length(x)-1))
  
  # loop through columns and calculate AUC
  for (i in 1:(length(x)-1)) {
    print(i)
    auc <- (x[i+1]-x[i])*(df_fraction[,i]) + (x[i+1]-x[i])*(df_fraction[,i+1]-df_fraction[,i])*0.5
    m_auc[,i] <- auc
  }
  
  # calculate AUC and add columns to data
  data_nomax$AUC <- round(rowSums(m_auc[,1:4])/4, 4)
  
  # make a min column specifying the minimum read per time point
  data_nomax$min <- rowMins(as.matrix(data_nomax[, c(paste0('N', 1:5))] + data_nomax[, c(paste0('F', 1:5))]))
  
  ## compare AUCs before and after removing the top flipping curve ----
  compareAUCs <- full_join(compareAUCs,
                           data_nomax %>%
                             select(TEVp, TEVs, AUC, min) %>%
                             rename(AUC_wo_max = AUC,
                                    min_wo_max = min),
                           by = c("TEVp", "TEVs"))
  
  # calculate the fold change in read-adjusted AUC that results from adding the
  # top plasmid curve to the other plasmid curves
  compareAUCs <- compareAUCs %>%
    # throw out samples where the AUC_wo_max is exactly 0 because this only happens
    # when all other replicates have too little reads to calculate an AUC
    filter(AUC_wo_max != 0) %>%
    mutate(fc = AUC / AUC_wo_max,
           # add a column specifying wheter the min without the top curve still reaches 20
           # (by definition, this is the case in the alternative min20DNA but not necessarily
           # for alternative min20protein)
           min_wo_max_over20 = min_wo_max >= 20) %>%
    arrange(desc(fc))
  
  # consider only cases where after removing the top curve the remaining curves
  # together still reach a minimum of 20 reads per time point (this step is important 
  # for the alternative min20protein)
  
  compareAUCs <- filter(compareAUCs, min_wo_max_over20 == TRUE)
  
  # save this table to a list
  list_compareAUCs[[alternative]] <- compareAUCs
  
  # # save this table to a file
  # write.table(compareAUCs, file = paste0("./02_DMS_controls_false_positives/results/compare_AUCs_", alternatives[alternative], ".txt"), sep = "\t", row.names=FALSE, quote = F)
  
  ### calculate how often false positives occur ----
  
  # initialize empty df
  df_false_freq <- data.frame()
  
  for (threshold in c(seq(1.1, 2.5, 0.1))) {
    frequency <- round(
      (compareAUCs %>%
         filter(fc >= threshold) %>% nrow()) /
        nrow(compareAUCs) * 100, 3)
    # write result in a table
    result <- data.frame(fc = threshold,
                         frequency = frequency)
    # bind this table with the initialized dataframe
    df_false_freq <- rbind(df_false_freq, result)
  }
  
  print(df_false_freq)
  
  # save result to file
  write.table(df_false_freq, file = paste0("./02_DMS_controls_false_positives/results/df_false_freq_", alternatives[alternative], ".txt"), sep = "\t", row.names=FALSE, quote = F)
  
  # create a density plot colored by TEVs
  p <- ggplot(compareAUCs, aes(x = fc, fill = TEVs)) +
    geom_density(alpha = 0.6) +
    annotate("text", x = 3, y = 10,
             label = paste("n =", num_unique_combinations)) +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "none") + # Remove the legend
    labs(fill = "TEVs") +
    ylab("Density") +
    coord_cartesian(ylim = c(0, 20)) #+
  #  scale_y_log10()  # Transform y-axis to a log scale
  
  # show plot
  print(p)
  # save plot
  save_plot(p, filename = paste0("./02_DMS_controls_false_positives/plots/", alternatives[alternative], "/false_positives_density.pdf"),
            base_height = 3,
            base_asp = 2)
  
}

## Plot the top cases of false positives ----------------------------------------

# Loop through the two alternative data sets
for (alternative in c(1:2)) {
  
  # get the correct data sets for the respective alternative
  compareAUCs <- list_compareAUCs[[alternative]]
  data <- list_data[[alternative]]
  data_BCs <- list_data_BCs[[alternative]]
  biol_reps <- list_biol_reps[[alternative]]
  
  for (i in c(1:5)) {
    # get TEVp variant
    hitTEVp <- compareAUCs$TEVp[i]
    # get TEVs variant
    hitTEVs <- compareAUCs$TEVs[i]
    # generate plot data
    # generate plot data
    p_TEVs_data <- data %>%
      # filtering
      filter(min >= 20,
             TEVp == hitTEVp) %>%
      # selecting only necessary columns
      select(TEVp, TEVs, t1, t2, t3, t4, t5) %>%
      # rename for correct plotting
      rename("0" = "t1") %>% rename("1" = "t2") %>% rename("2" = "t3") %>%
      rename("3" = "t4") %>% rename("4" = "t5") %>%
      # reformatting for plotting
      pivot_longer(cols = !c(TEVp, TEVs) , names_to = "timepoint", values_to = "fraction_flipped")
    
    # specify which TEVs to color
    p_TEVs_data$color <- ifelse(p_TEVs_data$TEVs == hitTEVs, "hitTEVs",
                                ifelse(p_TEVs_data$TEVs == "ENLYFQ*", "ENLYFQ*",
                                       ifelse(p_TEVs_data$TEVs == "ENLYFQS", "ENLYFQS", "other")))
    
    # Define custom legend labels, this is necessary to display the underlying value of hitTEVs in the plot
    legend_labels <- c("ENLYFQ*" = "ENLYFQ*", "ENLYFQ*" = "ENLYFQ*", "other" = "other", "hitTEVs" = hitTEVs)
    
    # generate plot data
    p_TEVs <- ggplot(p_TEVs_data, aes(x = as.integer(timepoint), y = fraction_flipped, col = color)) +
      geom_point() +
      geom_line(aes(group = TEVs)) +
      theme_cowplot(font_size = 12) +
      scale_x_continuous("Time after induction (h)", limits = c(0, 4)) +
      scale_y_continuous("Fraction flipped") +
      coord_cartesian(ylim = c(0, 0.6)) +
      labs(col = "TEVs") +
      ggtitle(paste("TEVp 0", hitTEVp)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = c("ENLYFQ*" = "#404040", "ENLYFQS" = "#0571b0", "hitTEVs" = "#ca0020", "other" = "#bababa"),
                         labels = legend_labels) +
      # additionally plot the hitTEVs curve again to ensure it is in the front
      geom_point(data = filter(p_TEVs_data, TEVs == hitTEVs), aes(group = TEVs), color = "#ca0020", size = 2) +
      geom_line(data = filter(p_TEVs_data, TEVs == hitTEVs), aes(group = TEVs), color = "#ca0020")
    
    
    ### go to the BC data set and extract this hit with its TEVs
    
    # First, calculate the fractions flipped for all BCs combined to add it to the plot later
    BCs_combined <- data_BCs %>%
      filter(TEVp == hitTEVp, TEVs == hitTEVs) %>%
      group_by(TEVp) %>%
      summarize(N1 = sum(N1), F1 = sum(F1),
                N2 = sum(N2), F2 = sum(F2),
                N3 = sum(N3), F3 = sum(F3),
                N4 = sum(N4), F4 = sum(F4),
                N5 = sum(N5), F5 = sum(F5)) %>%
      ungroup() %>%
      mutate(plasmidID = "combined",
             t1 = F1/sum(N1, F1),
             t2 = F2/sum(N2, F2),
             t3 = F3/sum(N3, F3),
             t4 = F4/sum(N4, F4),
             t5 = F5/sum(N5, F5),
             sum = sum(N1, N2, N3, N4, N5, F1, F2, F3, F4, F5),
             style = '5') %>% # add a style for plotting
      select(plasmidID, t1, t2, t3, t4, t5, sum, style)
    
    # prepare data for plotting
    hit <- data_BCs %>%
      # filtering
      filter(TEVp == hitTEVp, TEVs == hitTEVs) %>%
      mutate(plasmidID = paste0(BC_TEVp, "_", BC_TEVs),
             style = '1') %>% # add a style for plotting
      # selecting only necessary columns
      select(plasmidID, t1, t2, t3, t4, t5, sum, style) %>%
      # add fractions flipped for all BCs combined
      rbind(., BCs_combined) %>%
      # rename for correct plotting
      rename("0" = "t1") %>% rename("1" = "t2") %>% rename("2" = "t3") %>%
      rename("3" = "t4") %>% rename("4" = "t5") %>%
      # reformat for plotting
      pivot_longer(cols = c("0", "1", "2", "3", "4"), names_to = "timepoint", values_to = "fraction_flipped")
    
    
    # plot the TEVpBC plot
    p_TEVp <- ggplot(hit, aes(x = as.integer(timepoint), y = fraction_flipped, color = sum, group = plasmidID, linetype = style)) +
      theme_cowplot(font_size = 12) +
      geom_point() +
      geom_line() +
      scale_x_continuous("Time after induction (h)", limits = c(0, 4)) +
      scale_y_continuous("Fraction flipped") +
      coord_cartesian(ylim = c(0, 0.6)) +
      labs(col = "Total \nreads", linetype = "") + # change legend titles
      scale_color_gradient(low="blue", high="red") + # change color scale
      scale_linetype_discrete(labels = c("Individual \nplasmids", "Read \nadjusted \naverage")) + # change legend labels
      ggtitle(paste0("TEVs: ", hitTEVs)) + # add plot title
      theme(plot.title = element_text(hjust = 0.5)) +
      annotate("text", x = 3, y = 0.6,
               label = paste("fc =",
                             compareAUCs %>% filter(TEVp == hitTEVp, TEVs == hitTEVs) %>% pull(fc) %>% round(2)))
    
    # put the two plots next to each other
    pg <- plot_grid(p_TEVs, p_TEVp, labels = NULL)
    # show plot grid
    print(pg)
    # and save
    save_plot(pg, filename = paste0("./02_DMS_controls_false_positives/plots/", alternatives[alternative], "/fc_top", i, "_", hitTEVp, "_", hitTEVs, ".pdf"),
              base_height = 3,
              base_asp = 2.7)
  }
}


## Correct false positives ----

# Do the following again for both alternatives...

# define a cutoff for the fold change
fc_cutoff <- 2

# initialize a list to store the false positives in
list_false_positives <- list()

# initialize a list to store the results for replicates in
list_over3BCs_fc <- list()

# initialize a list to store the final results in
list_data_cleaned_by_fc <- list()

for (alternative in c(1:2)) {
  
  # get the correct data sets for the respective alternative
  compareAUCs <- list_compareAUCs[[alternative]]
  data <- list_data[[alternative]]
  data_BCs <- list_data_BCs[[alternative]]
  biol_reps <- list_biol_reps[[alternative]]
  
  # define false positives, i.e. plasmid IDs causing an fc above the cutoff
  falsepos <- compareAUCs %>%
    filter(fc >= fc_cutoff) %>%
    select(TEVp, TEVs, plasmidID_AUCmax) %>%
    rename(plasmidID = plasmidID_AUCmax) %>%
    mutate(TEVp_TEVs = paste0(TEVp, "_", TEVs))
  
  # save the false positives in a list
  list_false_positives[[alternative]] <- falsepos
  
  cat(paste0("Number of false positive plasmids (TEVp BC - TEVs BC combination) kicked out checking replicates using the alternative ", alternatives[[alternative]],":"),
      nrow(falsepos))
  
  # get the data of replicates, i.e. the underlying BC data
  biol_reps <- list_biol_reps[[alternative]]
  
  # filter out the false positive plasmids from the BC data
  biol_reps <- biol_reps %>% filter(!plasmidID %in% falsepos$plasmidID)
  
  #### recalculate AUCs without false positives ----
  
  # merge reads with the same TEVs_trans
  over3BCs_fc <- biol_reps %>%
    group_by(TEVp, TEVs) %>%
    summarize(
      replicates = n(),
      N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
      F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
      sum = sum(sum)) %>%
    ungroup()
  
  # add min
  over3BCs_fc$min <- rowMins(as.matrix(over3BCs_fc[, c(paste0('N', 1:5))] + over3BCs_fc[, c(paste0('F', 1:5))]))
  
  # filter
  over3BCs_fc <- over3BCs_fc %>% filter(min >= 20)
  
  # calculate fraction flipped
  for (tp in c(1:5)) {
    over3BCs_fc[, paste0('t', tp)] <- round(over3BCs_fc[, paste0('F', tp)] / 
                                              (over3BCs_fc[, paste0('F', tp)] + over3BCs_fc[, paste0('N', tp)]), 4)
  }
  
  # calculate AUC
  tps <- c(0, 60, 120, 180, 240)
  
  # calculate AUC
  x <- as.character(tps) # get time points
  df_y <- over3BCs_fc[, paste0('t', 1:5)] %>% as.data.frame(.) # get y values
  x <- as.integer(x)
  
  # initialize helper matrix for AUC calculation
  m_auc <- matrix(nrow = nrow(df_y), ncol = (length(x)-1))
  
  # loop through columns and calculate AUC
  for (i in 1:(length(x)-1)) {
    print(i)
    auc <- (x[i+1]-x[i])*(df_y[,i]) + (x[i+1]-x[i])*(df_y[,i+1]-df_y[,i])*0.5
    m_auc[, i] <- auc
  }
  
  over3BCs_fc$AUC <- round(rowSums(m_auc / 240), 4)
  
  # add important columns (n_mut_TEVs)
  over3BCs_fc <- inner_join(over3BCs_fc,
                            data_BCs %>% select(TEVp, TEVs, n_mut_TEVs) %>% distinct(),
                            by = c("TEVp", "TEVs"))
  
  # add information if the TEVp-TEVs combination was corrected by the fc
  over3BCs_fc <- over3BCs_fc %>%
    mutate(TEVp_TEVs = paste0(TEVp, "_", TEVs)) %>%
    mutate(source = ifelse(TEVp_TEVs %in% falsepos$TEVp_TEVs, "repl_fc_corrected", "repl_uncorrected")) %>%
    select(-TEVp_TEVs)
  
  # add the original AUC (i.e. read adjusted avg)
  over3BCs_fc <- left_join(over3BCs_fc,
                           data %>% select(TEVp, TEVs, AUC),
                           by = c("TEVp", "TEVs"),
                           suffix = c("", "_orig")) %>%
    arrange(desc(AUC)) %>%
    mutate(relative_change = AUC/AUC_orig)
  
  # # write table to file
  # write.table(over3BCs_fc, file = paste0("./02_DMS_controls_false_positives/results/replicates_cleaned_by_fc_", alternatives[[alternative]],".txt"), 
  #             sep = "\t", row.names=FALSE, quote = F)
  # save to list
  list_over3BCs_fc[[alternative]] <- over3BCs_fc
  
  # count unique TEVp-TEVs combinations
  replicates_fc_cleaned <- over3BCs_fc %>%
    # filter out any double mutants
    filter(!grepl("_", TEVp),
           !grepl("_", TEVs)) %>%
    nrow()
  cat(paste0("Number of unique TEVp-TEVs combinations (single mutants) with >= 3 replicates (", alternatives[[alternative]], ") after AUC correction by fc:"),
      replicates_fc_cleaned)
  
  # count unique TEVp-TEVs combinations initially
  replicates_fc <- biol_reps %>%
    # filter out any double mutants
    filter(!grepl("_", TEVp),
           !grepl("_", TEVs)) %>%
    mutate(TEVp_TEVs = paste0(TEVp, "_", TEVs)) %>%
    pull(TEVp_TEVs) %>%
    unique() %>%
    length()
  cat(paste0("Number of unique TEVp-TEVs combinations (single mutants) with >= 3 replicates (", alternatives[[alternative]], ") after AUC correction by fc:"),
      replicates_fc)    
  
  # calculate ratio of numbers
  ratio <- replicates_fc_cleaned / replicates_fc
  
  cat(paste0("Number of unique TEVp-TEVs combinations (single mutants) with >= 3 replicates (", alternatives[[alternative]], ") after vs before AUC correction by fc:"),
      ratio) 
  
  
}


# Generate final data sets   --------------------------------------------------

# I.e. without diverging duplicates and with corrected false positives


for (alternative in c(1:2)) {
  
  # get the correct data sets for the respective alternative
  data <- list_data[[alternative]]
  over3BCs_fc <- list_over3BCs_fc[[alternative]]
  diverging_duplicates <- list_diverging_duplicates[[alternative]]
  
  ## put together the corrected replicates and the rest of samples
  
  # get all other samples do add to the corrected replicates
  data_to_add <- anti_join(data, 
                           over3BCs_fc, by = c("TEVp", "TEVs")) %>%
    mutate(source = "original")
  
  # append filtered original samples to fc-corrected replicate samples
  data_cleaned_by_fc <- rbind(over3BCs_fc %>% select(names(data_to_add)),
                              data_to_add)
  
  ## remove diverging duplicates identified before
  data_cleaned_by_fc <- data_cleaned_by_fc %>%
    mutate(TEVp_TEVs = paste0(TEVp, "_", TEVs)) %>%
    filter(!TEVp_TEVs %in% diverging_duplicates)
  
  ## Lastly, throw out everything above the pos. ctrl.
  
  # get positive control
  pos_ctrl <- data_cleaned_by_fc %>%
    filter(TEVp == "0 C151A no ssrA",
           TEVs == "ENLYFQ*") %>%
    pull(AUC)
  
  # identify samples that yield a higher AUC than the positive control
  above_posctrl <- data_cleaned_by_fc %>%
    # ignore substrates containing a stop codon
    filter(!grepl("\\*", TEVs)) %>%
    filter(AUC >= pos_ctrl)
  
  # throw out these samples from the data
  data_cleaned_by_fc <- anti_join(data_cleaned_by_fc, above_posctrl)
  
  # write final data to file
  write.table(data_cleaned_by_fc, file = paste0("./02_DMS_controls_false_positives/results/Data_DMS_controls_fc_corrected_", alternatives[[alternative]],".txt"), 
              sep = "\t", row.names=FALSE, quote = F)
  
  # count unique TEVp-TEVs combinations of final data frame
  number_final_fc <- data_cleaned_by_fc %>%
    # filter out any double mutants
    filter(!grepl("_", TEVp),
           !grepl("_", TEVs)) %>%
    nrow
  cat(paste0("Total number of unique TEVp-TEVs combinations (single mutants) after AUC correction of replicates (", alternatives[[alternative]], ") by fc:"),
      number_final_fc)
  
  # save to list
  list_data_cleaned_by_fc[[alternative]] <- data_cleaned_by_fc
}

# Histogram of replicates after outlier correction------------------------------

for (alternative in c(1:2)) {
  
  # get the correct data sets for the respective alternative
  data_BCs <- list_data_BCs[[alternative]]
  data_cleaned_by_fc <- list_data_cleaned_by_fc[[alternative]]
  
  # count the number of replicates per TEVp-TEVs combination
  replicate_samples <- data_BCs %>%
    # for this plot, samples should be considered replicates only if their coverage
    # was high enough to calculate an AUC
    filter(is.finite(AUC)) %>%
    mutate(TEVp_TEVs = paste0(TEVp, "_", TEVs)) %>%
    filter(TEVp_TEVs %in% data_cleaned_by_fc$TEVp_TEVs) %>%
    group_by(TEVp, TEVs, TEVp_TEVs) %>% 
    summarize(n_replicates = n()) %>%
    # filter out wt TEVp and wt TEVs because they are extremely overrepresented
    filter(TEVp != "wt",
           TEVs != "ENLYFQS") %>%
    arrange(desc(n_replicates))
  
  # make a histogram plot of replicates
  p <- ggplot(replicate_samples, aes(x = n_replicates)) +
    geom_histogram(binwidth = 1, fill = "#8da0cb", color = "black", alpha = 0.7) +
    scale_x_continuous(breaks = unique(replicate_samples$n_replicates)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme_cowplot(font_size = 12) +
    annotate("text", x = 8, y = Inf, vjust = 1,
             label = paste("Average no. of replicates \nper sample =", round(sum(replicate_samples$n_replicates)/length(replicate_samples$TEVp_TEVs), 1))) +
    labs(x = "Number of replicates", y = "Count of unique TEVp-TEVs combinations") +
    geom_text(stat = 'count', aes(label = ..count..), size = 3, angle = 90, vjust = 0.5, hjust = -0.2)
  # print plot  
  print(p)
  
  #save plot
  save_plot(p, filename = paste0("./02_DMS_controls_false_positives/plots/histogram_replicates_fc_corrected_", alternatives[[alternative]], ".pdf"),
            base_height = 3.71, 
            base_asp = 1)
}



# ggplot(replicate_samples, aes(x = n_replicates)) +
#   geom_histogram(binwidth = 1, fill = "#8da0cb", color = "black", alpha = 0.7) +
#   scale_x_continuous(breaks = unique(replicate_samples$n_replicates)) +
#   scale_y_continuous(trans = 'log2') +
#   theme_cowplot(font_size = 12) +
#   labs(x = "Number of replicates", y = "Count of unique TEVp-TEVs combinations") +
#   geom_text(stat = 'count', aes(label = ..count..), size = 3, vjust = -0.5, hjust = 0.3)

for (alternative in c(1:2)) {
  cat(paste("Number of unique TEVp-TEVs combinations using alternative", alternatives[[alternative]], ":",
            list_data_cleaned_by_fc[[alternative]] %>% pull(TEVp_TEVs) %>% unique() %>% length(), "\n"))
}



# Histogram of BC distributions -------------------------------------------

for (alternative in c(1:2)) {
  
  # get the correct data sets for the respective alternative
  data_BCs <- list_data_BCs[[alternative]]
  diverging_duplicates <- list_diverging_duplicates[[alternative]]
  falsepos <- list_false_positives[[alternative]]
  
  # get the BC data after duplicate and outlier correction
  data_BCs_cleaned <- data_BCs %>%
    mutate(TEVp_TEVs = paste0(TEVp, "_", TEVs),
           plasmidID = paste0(BC_TEVp, "_", BC_TEVs)) %>%
    # remove the duplicates that diverge greatly from each other
    filter(!TEVp_TEVs %in% diverging_duplicates) %>%
    # remove the false positives
    filter(!plasmidID %in% falsepos)
  
  # summarize TEVp BC counts
  TEVp_BCs_counts <- data_BCs_cleaned %>%
    group_by(TEVp) %>% 
    summarize(n = n_distinct(BC_TEVp)) %>%
    # filter out wt TEVp
    filter(TEVp != "wt")
  
  # summarize TEVs BC counts
  TEVs_BCs_counts <- data_BCs_cleaned %>%
    group_by(TEVs) %>% 
    summarize(n = n_distinct(BC_TEVs)) %>%
    # filter out wt TEVs
    filter(TEVs != "ENLYFQS")
  
  # make a histogram plot of TEVs BC counts
  p <- ggplot(TEVs_BCs_counts, aes(x = n)) +
    geom_histogram(binwidth = 1, fill = "#bc9669", color = "black", alpha = 1) +
    # scale_y_continuous(breaks = seq(0, 35, by = 5), expand = c(0, 0)) +
    theme_cowplot(font_size = 12) +
    labs(x = "Number of TEVs barcodes", y = "Count of TEVs variants") +
    annotate("text", x = 10, y = Inf, vjust = 1,
             label = paste("Average no. of BCs \nper TEVs =", round(sum(TEVs_BCs_counts$n)/length(TEVs_BCs_counts$TEVs), 1))) 
  # print plot  
  print(p)
  #save plot
  save_plot(p, filename = paste0("./02_DMS_controls_false_positives/plots/histogram_TEVs_BCs_fc_corrected_", alternatives[[alternative]], ".pdf"),
            base_height = 3.71, 
            base_asp = 1)
}

