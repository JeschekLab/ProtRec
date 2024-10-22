# load libraries
library(tidyverse)
library(viridis)
library(cowplot)
library(stringi)
library(matrixStats)


# create folders ----------------------------------------------------------

# create a folder to store the individual plots and data in
dir.create("./01.1_DMS_mutated_backbones")
dir.create("./01.1_DMS_mutated_backbones/results")
dir.create("./01.1_DMS_mutated_backbones/plots")


# prepare data ------------------------------------------------------------

# read BC data before mapping
data_BCs <- read.table("./data/Data_DMS_unmapped.txt", sep = "\t", header = TRUE)

# throw out columns not needed
data_BCs <- data_BCs %>%
  select(-c(flask, indel, list_TEVp, reason_TEVp, source_TEVp, list_TEVs, reason_TEVs, source_TEVs))

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

data_BCs$AUC <- round(rowSums(m_auc[,1:4])/4, 4)


# Mutated TEVs backbones --------------------------------------------------


# Search for TEVs BCs (i.e. individual plasmid backbones) that consistently yield
# high AUCs irrespective of the TEVp variant

## inspection of TEVs BCs ----

# filter for min 20 reads and exclude the wt TEVs BC, i.e. the plasmid with wt TEVs
# that was used to clone the TEVp library
my_data <- data_BCs%>%
  filter(min >= 20) %>%
  select(BC_TEVp, BC_TEVs, AUC) %>%
  # throw out wt TEVs BC
  filter(BC_TEVs != "ACGTCGCTGA")

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
            aes(x = factor(BC_TEVs, levels=TEVs_BC_order), y = AUC, fill = n, group = BC_TEVs)) +
  geom_boxplot(width = 0.8, color = "black", alpha = 1) +
  labs(x = "TEVs BC (= plasmid backbone ID)", y = "AUC") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# print plot
print(p)
#save plot
save_plot(p, filename = "./01.1_DMS_mutated_backbones/plots/TEVs_BC_backbones_highest.pdf",
          base_height = 3.71, 
          base_asp = 2)

# Aligning the raw Pacbio reads of the 20 TEVs BCs with the highest median AUC
# revealed that the first 11 TEVs BCs showed mutations
# list the mutated TEVs BCs in a table
mutated_backbones <- data_frame(BC_TEVs = TEVs_BC_medians$BC_TEVs[1:11],
                                identified_by = "high AUC median")

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
            aes(x = factor(BC_TEVs, levels=TEVs_BC_order), y = AUC, fill = n,
                group = BC_TEVs)) +
  geom_boxplot(width = 0.8, color = "black", alpha = 1) +
  labs(x = "TEVs BC (= plasmid backbone ID)", y = "AUC") +
  theme_cowplot(font_size = 12) +
  coord_cartesian(ylim = c(0, 0.25)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# print plot
print(p)
#save plot
save_plot(p, filename = "./01.1_DMS_mutated_backbones/plots/TEVs_BC_backbones_lowest.pdf",
          base_height = 3.71, 
          base_asp = 2)


# Aligning the raw Pacbio reads of some TEVs BCs with the lowest median AUC
# revealed that at least the first 30 TEVs BCs showed mutations --> 
# I have no clue how many more I would need to inspect by manual alignment -->
# I will plot and inspect all TEVs variants individually:


## inspection of TEVs variants separately ----

# overview of how many BCs we have per TEVs
TEVs_summary <- data_BCs %>% 
  filter(n_mut_TEVs <= 1) %>%
  select(BC_TEVs, TEVs) %>%
  distinct() %>%
  group_by(TEVs) %>%
  summarize(n_TEVs_BCs = n()) %>%
  arrange(desc(n_TEVs_BCs))

# create folder to store plots in
dir.create("./01.1_DMS_mutated_backbones/plots/TEVs_BCs_individual")

# # make a box plot for each TEVs
# for (i in TEVs_summary$TEVs) {
#   p <- ggplot(data = my_data %>% filter(TEVs == i),
#          aes(x = reorder(BC_TEVs, AUC, median), y = AUC, fill = n)) +
#     # geom_violin() +
#     geom_boxplot(width = 0.8, color = "black", alpha = 1) +
#     labs(x = "TEVs BC (= plasmid backbone ID)", y = "AUC") +
#     theme_cowplot(font_size = 12) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),
#           plot.title = element_text(hjust = 0.5)) +
#     ggtitle(paste("TEVs:", i)) # add plot title
#   print(p)
# 
#   #save plot
#   save_plot(p, filename = paste0("./01.1_DMS_mutated_backbones/plots/TEVs_BCs_individual/TEVs_BCs_", i, ".pdf"),
#             base_height = 3.71,
#             base_asp = 2)
#   }

# write down suspicious TEVs BCs in a table
# read the table of suspicious TEVs BCs identified through a visual inspection of boxplots
TEVs_BC_manual <- read.table("./data/TEVs_BCs_suspicious.txt", sep = "\t", header = FALSE) %>%
  rename(BC_TEVs = V1)
# clean the data: remove white space from manually copy pasted TEVs BCs
TEVs_BC_manual$BC_TEVs <- trimws(TEVs_BC_manual$BC_TEVs)

# add TEVs column
TEVs_BC_manual <- left_join(TEVs_BC_manual,
                            data_BCs %>% select(BC_TEVs, TEVs) %>% distinct())

# Aligning the raw Pacbio reads of the manually selected suspicious TEVs BCs
# confirmed a few TEVs BCs with mutations --> kick them out
confirmed_correct <- c("AAGTCGCTGA", "ACGTCGCTTA", "ACTTCGCTGA", "CCCCCCCCTG",
                       "GAAAAGGTTA", "GCCCCTCCCC", "GGGGGTGGC", "CCGCCCCACG",
                       "TCGTGGGATC")

# define mutated TEVs BCs yielding consitently low AUCs
mutatedTEVsBCs <- TEVs_BC_manual %>% filter(!BC_TEVs %in% confirmed_correct) %>% pull(BC_TEVs)
# add one more mutated backbone that I found by coincident when investigating the controls data set
mutatedTEVsBCs <- c(mutatedTEVsBCs, "GAATAAGAAA")

# add these mutated TEVs BCs to the table
mutated_backbones <- rbind(mutated_backbones,
                           data_frame(BC_TEVs = mutatedTEVsBCs,
                                      identified_by = "individual inspection of TEVs BC box plots"))

# save the table of mutated TEVs BCs (i.e., plasmid backbones)
write.table(mutated_backbones, file = "./01.1_DMS_mutated_backbones/results/mutated_backbones.txt", sep = "\t", row.names=FALSE, quote = F)



# Save data ----------------------------------------------------


# throw out all identified mutated TEVs BCs from the data
data_BCs <- filter(data_BCs, !BC_TEVs %in% mutated_backbones$BC_TEVs)

# save the curated data of unmapped TEVs BCs
write.table(data_BCs, file = "./01.1_DMS_mutated_backbones/results/Data_DMS_unmapped_without_mutated_backbones.txt", sep = "\t", row.names=FALSE, quote = F)


