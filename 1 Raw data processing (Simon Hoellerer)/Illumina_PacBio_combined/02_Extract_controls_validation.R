# import variables
ROOT_DIR <- '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/uASPIre_screen4_5/Illumina_PacBio_combined'
PLT_DIR <- paste0(ROOT_DIR, '/plots')
RES_DIR <- paste0(ROOT_DIR, '/results')
Illumina_DIR <- '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/uASPIre_screen4_5/2022_LH_SH_S4_combined'
PacBio_DIR <- '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio_combined'

# change directory
setwd(ROOT_DIR)

# load libraries
library(tidyverse)
library(VennDiagram)
library(spgs)
library(data.table)
library(stringdist)
library(matrixStats)
library(ggseqlogo)
library(bioseq)

# create directories
dir.create(RES_DIR, showWarnings = FALSE)
dir.create(PLT_DIR, showWarnings = FALSE)

# create color palette
COLORS <- c('#F8766D', '#D89000', '#A3A500', '#39B600', '#00BF7D', '#00BFC4',
  '#00B0F6', '#9590FF', '#E76BF3', '#FF62BC')
COLUMNS <- c(paste0('N', 1:5), paste0('F', 1:5))

# load Illumina data
raw_DMS <- read.table(
  paste0(Illumina_DIR, '/results/data_DMS.txt'),
  header = T, sep = '\t',
  colClasses = c(rep('character', 2), rep('integer', 10), 'character')) %>%
  rename(flask = library)

raw_MLD <- read.table(
  paste0(Illumina_DIR, '/results/data_MLD.txt'),
  header = T, sep = '\t',
  colClasses = c(rep('character', 2), rep('integer', 10), 'character')) %>%
  rename(flask = library)

### here ###

# load PacBio data
ref_TEVp_MLD <- read.table(
  paste0(PacBio_DIR, '/TEVp/results/BC_TEVp_PacBio_MLD_clean.txt'),
  header = T, sep = '\t',
  colClasses = 'character') %>%
  rename(list_TEVp = list) %>%
  rename(reason_TEVp = reason) %>%
  mutate(source_TEVp = 'PacBio')

ref_TEVp_DMS <- read.table(
  paste0(PacBio_DIR, '/TEVp/results/BC_TEVp_PacBio_DMS_clean.txt'),
  header = T, sep = '\t',
  colClasses = 'character') %>%
  rename(list_TEVp = list) %>%
  rename(reason_TEVp = reason) %>%
  mutate(source_TEVp = 'PacBio') %>%
  rename(n_mut_TEVp = n_mut)

ref_TEVs_DMS <- read.table(
  paste0(PacBio_DIR, '/TEVs/results/BC_TEVs_PacBio_DMS_clean.txt'),
  header = T, sep = '\t',
  colClasses = 'character') %>%
  rename(list_TEVs = list) %>%
  rename(reason_TEVs = reason) %>%
  mutate(source_TEVs = 'PacBio') %>%
  rename(n_mut_TEVs = n_mut)

# load manual cloned and known barcodes
ref_TEVs_controls <- read.table(
  paste0(ROOT_DIR, '/BC_TEVs_controls.txt'),
  header = T, sep = '\t') %>%
  mutate(source_TEVs = 'controls')

# load manual cloned and known barcodes
ref_TEVs_validation <- read.table(
  paste0(ROOT_DIR, '/BC_TEVs_validation.txt'),
  header = T, sep = '\t') %>%
  mutate(source_TEVs = 'validation')

ref_TEVp_controls <- read.table(
  paste0(ROOT_DIR, '/BC_TEVp_controls.txt'),
  header = T, sep = '\t') %>%
  mutate(source_TEVp = 'controls')

ref_TEVp_validation <- read.table(
  paste0(ROOT_DIR, '/BC_TEVp_validation.txt'),
  header = T, sep = '\t') %>%
  mutate(source_TEVp = 'validations')

ref_NNK_theoretical <- read.table(
  paste0(ROOT_DIR, '/NNK_theoretical.txt'),
  header = T, sep = '\t') %>%
  mutate(AA = seq_translate(dna(DNA)))

# add sum to data
raw_MLD$sum <- rowSums(raw_MLD[, COLUMNS])
raw_DMS$sum <- rowSums(raw_DMS[, COLUMNS])

##################### DMS INTERNAL CONTROLS #####################
# extract library
data <- raw_DMS %>%
  inner_join(., ref_TEVp_controls, by = 'BC_TEVp') %>%
  inner_join(., ref_TEVs_controls, by = 'BC_TEVs')

# merge reads with the same TEVs_trans (doesnt do anything here)
data <- data %>%
  group_by(TEVp, TEVs, PlasmidID, libraryID, Note, source_TEVp,
    source_TEVs) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add min
data$min <- rowMins(as.matrix(data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))]))

# filter
data <- data %>% filter(min >= 20)

# calculate fraction flipped
for (tp in c(1:5)) {
  data[, paste0('t', tp)] <- round(data[, paste0('F', tp)] / 
    (data[, paste0('F', tp)] + data[, paste0('N', tp)]), 4)
}

# calculate AUC
tps <- c(0, 60, 120, 180, 240)

# calculate AUC
x <- as.character(tps) # get time points
df_y <- data[, paste0('t', 1:5)] %>% as.data.frame(.) # get y values
x <- as.integer(x)

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_y), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_y[,i]) + (x[i+1]-x[i])*(df_y[,i+1]-df_y[,i])*0.5
  m_auc[, i] <- auc
}

data$AUC <- round(rowSums(m_auc / 240), 4)

# write to file
write.table(data,
  file = paste0(RES_DIR, '/raw_DMS_internal_controls.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')



##################### DMS CONTROLS #####################
# extract library
data <- raw_DMS %>%
  inner_join(., ref_TEVp_controls, by = 'BC_TEVp') %>%
  inner_join(., ref_TEVs_DMS %>% filter(list_TEVs == 'WL'), by = 'BC_TEVs')

# merge reads with the same TEVs_trans (doesnt do anything here)
data <- data %>%
  group_by(TEVp, TEVs, n_mut_TEVs, PlasmidID, libraryID, Note, source_TEVp,
    source_TEVs) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add min
data$min <- rowMins(as.matrix(data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))]))

# filter
data <- data %>% filter(min >= 20)

# calculate fraction flipped
for (tp in c(1:5)) {
  data[, paste0('t', tp)] <- round(data[, paste0('F', tp)] / 
    (data[, paste0('F', tp)] + data[, paste0('N', tp)]), 4)
}

# calculate AUC
tps <- c(0, 60, 120, 180, 240)

# calculate AUC
x <- as.character(tps) # get time points
df_y <- data[, paste0('t', 1:5)] %>% as.data.frame(.) # get y values
x <- as.integer(x)

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_y), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_y[,i]) + (x[i+1]-x[i])*(df_y[,i+1]-df_y[,i])*0.5
  m_auc[, i] <- auc
}

data$AUC <- round(rowSums(m_auc / 240), 4)

# write to file
write.table(data,
  file = paste0(RES_DIR, '/raw_DMS_controls.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')



##################### DMS VALIDATION MANUAL #####################
# extract library
data <- raw_DMS %>%
  inner_join(., ref_TEVp_validation, by = 'BC_TEVp') %>%
  inner_join(., ref_TEVs_validation, by = 'BC_TEVs')

# merge reads with the same TEVs_trans (doesnt do anything here)
data <- data %>%
  group_by(TEVp, TEVs, PlasmidID, libraryID, source_TEVp,
    source_TEVs) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add min
data$min <- rowMins(as.matrix(data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))]))

# filter
data <- data %>% filter(min >= 20)

# calculate fraction flipped
for (tp in c(1:5)) {
  data[, paste0('t', tp)] <- round(data[, paste0('F', tp)] / 
    (data[, paste0('F', tp)] + data[, paste0('N', tp)]), 4)
}

# calculate AUC
tps <- c(0, 60, 120, 180, 240)

# calculate AUC
x <- as.character(tps) # get time points
df_y <- data[, paste0('t', 1:5)] %>% as.data.frame(.) # get y values
x <- as.integer(x)

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_y), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_y[,i]) + (x[i+1]-x[i])*(df_y[,i+1]-df_y[,i])*0.5
  m_auc[, i] <- auc
}

data$AUC <- round(rowSums(m_auc / 240), 4)

# write to file
write.table(data,
  file = paste0(RES_DIR, '/raw_DMS_validation_manual.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')



##################### DMS VALIDATION PACBIO #####################
# extract library
data <- raw_DMS %>%
  inner_join(., ref_TEVp_validation, by = 'BC_TEVp') %>%
  inner_join(., ref_TEVs_DMS %>% filter(list_TEVs == 'WL'), by = 'BC_TEVs')

# merge reads with the same TEVs_trans (doesnt do anything here)
data <- data %>%
  group_by(TEVp, TEVs, n_mut_TEVs, PlasmidID, libraryID, source_TEVp,
    source_TEVs) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add min
data$min <- rowMins(as.matrix(data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))]))

# filter
data <- data %>% filter(min >= 20)

# calculate fraction flipped
for (tp in c(1:5)) {
  data[, paste0('t', tp)] <- round(data[, paste0('F', tp)] / 
    (data[, paste0('F', tp)] + data[, paste0('N', tp)]), 4)
}

# calculate AUC
tps <- c(0, 60, 120, 180, 240)

# calculate AUC
x <- as.character(tps) # get time points
df_y <- data[, paste0('t', 1:5)] %>% as.data.frame(.) # get y values
x <- as.integer(x)

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_y), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_y[,i]) + (x[i+1]-x[i])*(df_y[,i+1]-df_y[,i])*0.5
  m_auc[, i] <- auc
}

data$AUC <- round(rowSums(m_auc / 240), 4)

# write to file
write.table(data,
  file = paste0(RES_DIR, '/raw_DMS_validation_PacBio.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')



##################### MLD INTERNAL CONTROLS #####################
# extract library
data <- raw_MLD %>%
  inner_join(., ref_TEVp_controls, by = 'BC_TEVp') %>%
  inner_join(., ref_TEVs_controls, by = 'BC_TEVs')

# merge reads with the same TEVs_trans (doesnt do anything here)
data <- data %>%
  group_by(TEVp, TEVs, PlasmidID, libraryID, Note, source_TEVp,
    source_TEVs) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add min
data$min <- rowMins(as.matrix(data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))]))

# filter
data <- data %>% filter(min >= 20)

# calculate fraction flipped
for (tp in c(1:5)) {
  data[, paste0('t', tp)] <- round(data[, paste0('F', tp)] / 
    (data[, paste0('F', tp)] + data[, paste0('N', tp)]), 4)
}

# calculate AUC
tps <- c(0, 60, 120, 180, 240)

# calculate AUC
x <- as.character(tps) # get time points
df_y <- data[, paste0('t', 1:5)] %>% as.data.frame(.) # get y values
x <- as.integer(x)

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_y), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_y[,i]) + (x[i+1]-x[i])*(df_y[,i+1]-df_y[,i])*0.5
  m_auc[, i] <- auc
}

data$AUC <- round(rowSums(m_auc / 240), 4)

# write to file
write.table(data,
  file = paste0(RES_DIR, '/raw_MLD_internal_controls.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')
