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

##################### plotting #####################
### plot overlap between PacBio and MLD ###

# merge MLD Illumina with PacBio
data <- raw_MLD %>%
  inner_join(., ref_TEVp_MLD %>% filter(list == 'WL'), by = 'BC_TEVp') %>%
  inner_join(., ref_TEVs_controls, by = 'BC_TEVs')

# sum reads for each TEVp_barcode
data_DNA <- data %>% group_by(BC_TEVp, library) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add rowmin and sort
data_DNA$min <- rowMins(as.matrix(data_DNA[, c(paste0('N', 1:5))] + data_DNA[, c(paste0('F', 1:5))]))
data_DNA <- data_DNA %>% arrange(desc(sum))

# filter by cutoff and write to file
MLD_number_of_TEVp_BC_variants <- data_DNA %>%
  filter(min >= 20) %>%
  group_by(library) %>%
  summarize(n = n())

write.table(MLD_number_of_TEVp_BC_variants,
  file = paste0(RES_DIR, '/MLD_number_of_TEVp_BC_variants.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

# sum reads for each TEVp
data_PRO <- data %>% group_by(TEVp_ABC, TEVp_DEF, TEVp_XYZ, library) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add rowmin and sort
data_PRO$min <- rowMins(as.matrix(data_PRO[, c(paste0('N', 1:5))] + data_PRO[, c(paste0('F', 1:5))]))
data_PRO <- data_PRO %>% arrange(desc(sum))

# filter by cutoff and write to file
MLD_number_of_TEVp_variants <- data_PRO %>%
  filter(min >= 20) %>%
  group_by(library) %>%
  summarize(n = n())

write.table(MLD_number_of_TEVp_variants,
  file = paste0(RES_DIR, '/MLD_number_of_TEVp_variants.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')



##################### plot overlap between PacBio and DMS - TEVp #####################
# merge DMS Illumina with PacBio
data <- raw_DMS %>%
  inner_join(., ref_TEVp_DMS %>% filter(list == 'WL'), by = 'BC_TEVp') %>%
  inner_join(., ref_TEVs_DMS %>% filter(list == 'WL'), by = 'BC_TEVs')

# select only single mutants
data <- data %>%
  filter(n_mut == 1)

# add mutation information
data <- data %>%
  mutate(codon = readr::parse_number(TEVp))

# filter
data <- data %>%
  filter(codon >= 1, codon <= 234)
  
# add bin information
data <- data %>%
  mutate(bin = ceiling(codon / 24))

# sum reads for each TEVp
data <- data %>% group_by(TEVp, bin, BC_TEVp) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add rowmin and sort
data$min <- rowMins(as.matrix(data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))]))
data <- data %>% arrange(desc(sum))

# count number of barcodes per variant
DMS_number_of_BC_per_TEVp <- data %>%
  group_by(TEVp, bin) %>%
  summarize(n = n()) %>%
  ungroup()

# make boxplots
p <- ggplot(DMS_number_of_BC_per_TEVp, aes(x = as.factor(bin), y = n, fill = as.factor(bin))) +
  stat_boxplot(geom = 'errorbar') +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.2, alpha = 0.5) +
  scale_x_discrete('bin') +
  scale_y_continuous('number of BCs per TEVp') +
  theme_SH() +
  theme(legend.position = 'none') +
  coord_cartesian(clip = 'off')

ggsave(plot = p, paste0(PLT_DIR, '/DMS_number_of_BC_per_TEVp.png'),
  width = 6, height = 3)

# write to file
write.table(DMS_number_of_BC_per_TEVp,
  file = paste0(RES_DIR, '/DMS_number_of_BC_per_TEVp.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

# sum reads for each TEVp
data_TEVp <- data %>% group_by(TEVp, bin) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add rowmin and sort
data_TEVp$min <- rowMins(as.matrix(data_TEVp[, c(paste0('N', 1:5))] + data_TEVp[, c(paste0('F', 1:5))]))
data_TEVp <- data_TEVp %>% arrange(desc(sum))

# filter by cutoff and write to file
DMS_number_of_TEVp_variants <- data_TEVp %>%
  filter(min >= 20) %>%
  group_by(bin) %>%
  summarize(n = n()) %>%
  ungroup()

write.table(DMS_number_of_TEVp_variants,
  file = paste0(RES_DIR, '/DMS_number_of_TEVp_variants.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

# sum reads for each TEVp_BC
data_TEVp_BC <- data %>% group_by(bin, BC_TEVp) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add rowmin and sort
data_TEVp_BC$min <- rowMins(as.matrix(data_TEVp_BC[, c(paste0('N', 1:5))] + data_TEVp_BC[, c(paste0('F', 1:5))]))
data_TEVp_BC <- data_TEVp_BC %>% arrange(desc(sum))

# filter by cutoff and write to file
DMS_number_of_TEVp_BC_variants <- data_TEVp_BC %>%
  filter(min >= 20) %>%
  group_by(bin) %>%
  summarize(n = n()) %>%
  ungroup()

write.table(DMS_number_of_TEVp_BC_variants,
  file = paste0(RES_DIR, '/DMS_number_of_TEVp_BC_variants.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')



##################### plot overlap between PacBio and DMS - TEVs #####################
# sum reads for each TEVs_barcode
data_DNA <- data %>% group_by(BC_TEVs) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add rowmin and sort
data_DNA$min <- rowMins(as.matrix(data_DNA[, c(paste0('N', 1:5))] + data_DNA[, c(paste0('F', 1:5))]))
data_DNA <- data_DNA %>% arrange(desc(sum))

# filter by cutoff and write to file
DMS_number_of_TEVs_BC_variants <- data_DNA %>%
  filter(min >= 20) %>%
  summarize(n = n())

DMS_number_of_TEVs_BC_variants <- write.table(DMS_number_of_TEVs_BC_variants,
  file = paste0(RES_DIR, '/DMS_number_of_TEVs_BC_variants.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

# sum reads for each TEVs
data_PRO <- data %>% group_by(TEVs) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add rowmin and sort
data_PRO$min <- rowMins(as.matrix(data_PRO[, c(paste0('N', 1:5))] + data_PRO[, c(paste0('F', 1:5))]))
data_PRO <- data_PRO %>% arrange(desc(sum))

# filter by cutoff
DMS_number_of_TEVs_variants <- data_PRO %>%
  filter(min >= 20) %>%
  summarize(n = n())

# write to file
write.table(DMS_number_of_TEVs_variants,
  file = paste0(RES_DIR, '/DMS_number_of_TEVs_variants.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')



##################### blacklist whitelist ratio - DMS TEVp #####################
### ratio of BL WL in DMS ###
data <- raw_DMS %>%
  left_join(., ref_TEVp_DMS, by = 'BC_TEVp') %>%
  inner_join(., ref_TEVs_DMS %>% filter(list == 'WL'), by = 'BC_TEVs')

sum(raw_DMS$sum) # 1978158932
sum(data$sum) # 1539258879

# combine reads by list
data <- data %>%
  group_by(BC_TEVp, list_TEVs, reason_TEVs) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add rowmin and sort
data$min <- rowMins(as.matrix(data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))]))
data <- data %>% arrange(desc(sum))

# filter by cutoff 
DMS_number_of_TEVp_BC_WLBL <- data %>%
  filter(min >= 20) %>%
  group_by(list_TEVs, reason_TEVs) %>%
  summarize(total_reads = sum(sum)) %>%
  ungroup() %>%
  filter(reason != 'synonymous')

names(DMS_number_of_TEVp_BC_WLBL) <- c('list', 'reason', 'total_reads')
DMS_number_of_TEVp_BC_WLBL$reason[DMS_number_of_TEVp_BC_WLBL$reason == 'none'] <- 'WL'
DMS_number_of_TEVp_BC_WLBL$reason[is.na(DMS_number_of_TEVp_BC_WLBL$reason)] <- 'only Illumina'

# plot
p <- ggplot(DMS_number_of_TEVp_BC_WLBL, aes(x = '', y = total_reads, fill = reason)) +
  geom_bar(stat = 'identity', width = 1) +
  coord_polar('y', start = 0) +
  theme_void() +
  scale_fill_manual(values = c(COLORS[c(9, 5, 10)]))

ggsave(paste0(PLT_DIR, '/DMS_number_of_TEVp_BC_WLBL.png'),
  plot = p, width = 3, height = 3)

# write to file
write.table(DMS_number_of_TEVp_BC_WLBL,
  file = paste0(RES_DIR, '/DMS_number_of_TEVp_BC_WLBL.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')



##################### blacklist whitelist ratio - DMS TEVs #####################
### ratio of BL WL in DMS ###
data <- raw_DMS %>%
  inner_join(., ref_TEVp_DMS %>% filter(list == 'WL'), by = 'BC_TEVp') %>%
  left_join(., ref_TEVs_DMS, by = 'BC_TEVs')

sum(raw_DMS$sum) # 1978158932
sum(data$sum)     # 1539258879

# combine reads by list
data <- data %>%
  group_by(BC_TEVs, list_TEVs, reason_TEVs) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add rowmin and sort
data$min <- rowMins(as.matrix(data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))]))
data <- data %>% arrange(desc(sum))

# filter by cutoff 
DMS_number_of_TEVs_BC_WLBL <- data %>%
  filter(min >= 20) %>%
  group_by(list_TEVs, reason_TEVs) %>%
  summarize(total_reads = sum(sum)) %>%
  ungroup()

names(DMS_number_of_TEVs_BC_WLBL) <- c('list', 'reason', 'total_reads')
DMS_number_of_TEVs_BC_WLBL$reason[DMS_number_of_TEVs_BC_WLBL$reason == 'none'] <- 'WL'
DMS_number_of_TEVs_BC_WLBL$reason[is.na(DMS_number_of_TEVs_BC_WLBL$reason)] <- 'only Illumina'

# plot
p <- ggplot(DMS_number_of_TEVs_BC_WLBL, aes(x = '', y = total_reads, fill = reason)) +
  geom_bar(stat = 'identity', width = 1) +
  coord_polar('y', start = 0) +
  theme_void() +
  scale_fill_manual(values = c(COLORS[c(9, 5, 10)]))

ggsave(paste0(PLT_DIR, '/TEVs_BL_WL.png'),
  plot = p, width = 3, height = 3)

# write to file
write.table(DMS_number_of_TEVs_BC_WLBL,
  file = paste0(RES_DIR, '/DMS_number_of_TEVs_BC_WLBL.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')



##################### blacklist whitelist ratio - MLD TEVp #####################
### ratio of BL WL in DMS ###
data <- raw_MLD %>%
  left_join(., ref_TEVp_MLD, by = 'BC_TEVp') %>%
  inner_join(., ref_TEVs_controls, by = 'BC_TEVs')

sum(raw_MLD$sum) # 2590265690
sum(data$sum)     # 2496014526

# combine reads by list
data <- data %>%
  group_by(BC_TEVs, list_TEVp, reason_TEVp) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
    sum = sum(sum)) %>%
  ungroup()

# add rowmin and sort
data$min <- rowMins(as.matrix(data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))]))
data <- data %>% arrange(desc(sum))

# filter by cutoff 
MLD_number_of_TEVp_BC_WLBL <- data %>%
  filter(min >= 20) %>%
  group_by(list_TEVp, reason_TEVp) %>%
  summarize(total_reads = sum(sum)) %>%
  ungroup()

MLD_number_of_TEVp_BC_WLBL <- MLD_number_of_TEVp_BC_WLBL %>%
  rename(list = list_TEVp, reason = reason_TEVp)

MLD_number_of_TEVp_BC_WLBL$reason[is.na(MLD_number_of_TEVp_BC_WLBL$reason)] <- 'only Illumina'

# plot
p <- ggplot(MLD_number_of_TEVp_BC_WLBL, aes(x = '', y = total_reads, fill = reason)) +
  geom_bar(stat = 'identity', width = 1) +
  coord_polar('y', start = 0) +
  theme_void() +
  scale_fill_manual(values = c(COLORS[c(5, 7, 8, 9, 10)]))

ggsave(plot = p, paste0(PLT_DIR, '/MLD_number_of_TEVp_BC_WLBL.png'),
  width = 3, height = 3)

# write to file
write.table(MLD_number_of_TEVp_BC_WLBL,
  file = paste0(RES_DIR, '/MLD_number_of_TEVp_BC_WLBL.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')



##################### expected distribution for ABC - 9k #####################
# directories
N9K_DIR <- '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/2021_LH_9k/results'
LUP_DIR <- '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/2021_LH_TEV_pop2/2021_LH_TEV_pop2_PCR/results'

# read lookup table
lookup <- read.table(paste0(LUP_DIR, '/02_LH_TEV_pop2_lookup_final_HQ.txt'),
  sep = '\t', header = TRUE, colClasses = 
  c(rep('character', 2), rep('integer', 3), 'numeric'))

# group reads
lookup <- lookup %>%
  filter(!(grepl('\\*', TEVp))) %>%
  separate(., col = TEVp,
    into = c('A29', 'A30', 'A31', 'A32'), sep = '', remove = T) %>%
  select(-A29) %>%
  group_by(A30, A31, A32) %>%
  summarize(sum = sum(reads_BC)) %>%
  ungroup()

A30 <- lookup %>% 
  group_by(A30) %>%
  summarize(sum30 = sum(sum)) %>%
  rename(AA = A30)

A31 <- lookup %>% 
  group_by(A31) %>%
  summarize(sum31 = sum(sum)) %>%
  rename(AA = A31)

A32 <- lookup %>% 
  group_by(A32) %>%
  summarize(sum32 = sum(sum)) %>%
  rename(AA = A32)

data_plot_lookup <- A30 %>%
  inner_join(., A31, by = 'AA') %>%
  inner_join(., A32, by = 'AA')

### same for MLD
lookup_MLD <- read.table(paste0(PacBio_DIR, '/TEVp/results/BC_TEVp_PacBio_MLD.txt'),
  sep = '\t', header = TRUE)

# filter and group
lookup_MLD <- lookup_MLD %>%
  filter(list == 'WL') %>%
  filter(library == 'ABCDEF') %>%
  group_by(TEVp_ABC) %>%
  summarize(sum = sum(reads_TEVpxBC)) %>%
  ungroup() %>%
  rename(TEVp = TEVp_ABC) %>%
  filter(!(grepl('\\*', TEVp))) %>%
  separate(., col = TEVp,
    into = c('A29', 'A30', 'A31', 'A32'), sep = '', remove = T) %>%
  select(-A29) %>%
  group_by(A30, A31, A32) %>%
  summarize(sum = sum(sum)) %>%
  ungroup()

A30 <- lookup_MLD %>% 
  group_by(A30) %>%
  summarize(sum30 = sum(sum)) %>%
  rename(AA = A30)

A31 <- lookup_MLD %>% 
  group_by(A31) %>%
  summarize(sum31 = sum(sum)) %>%
  rename(AA = A31)

A32 <- lookup_MLD %>% 
  group_by(A32) %>%
  summarize(sum32 = sum(sum)) %>%
  rename(AA = A32)

data_plot_lookup_MLD <- A30 %>%
  inner_join(., A31, by = 'AA') %>%
  inner_join(., A32, by = 'AA')

# correct for NNK bias
NNK <- ref_NNK_theoretical %>%
  group_by(AA) %>%
  summarize(n = n()) %>%
  ungroup()

data_plot_lookup <- data_plot_lookup %>%
  inner_join(., NNK, by = 'AA')
data_plot_lookup_MLD <- data_plot_lookup_MLD %>%
  inner_join(., NNK, by = 'AA')

# correct
data_plot_lookup <- data_plot_lookup %>%
  mutate(
    sum30 = sum30 / n,
    sum31 = sum31 / n,
    sum32 = sum32 / n) %>%
  mutate(sum30 = sum30 / sum(sum30),
    sum31 = sum31 / sum(sum31),
    sum32 = sum32 / sum(sum32)) %>%
  select(-n)

# correct
data_plot_lookup_MLD <- data_plot_lookup_MLD %>%
  mutate(
    sum30 = sum30 / n,
    sum31 = sum31 / n,
    sum32 = sum32 / n) %>%
  mutate(sum30 = sum30 / sum(sum30),
    sum31 = sum31 / sum(sum31),
    sum32 = sum32 / sum(sum32)) %>%
  select(-n)

ROWNAMES <- data_plot_lookup$AA
data_plot_lookup <- as.matrix(data_plot_lookup %>% select(-AA))
rownames(data_plot_lookup) <- ROWNAMES

ROWNAMES <- data_plot_lookup_MLD$AA
data_plot_lookup_MLD <- as.matrix(data_plot_lookup_MLD %>% select(-AA))
rownames(data_plot_lookup_MLD) <- ROWNAMES

p <- ggplot() +
  geom_logo(data = data_plot_lookup - 0.05, method = 'custom', seq_type = 'aa') +
  theme_SH() +
  scale_x_continuous('randomised position in TEVp',
    breaks = c(1, 2, 3), labels = 30:32) +
  scale_y_continuous('AA fraction',
    limits = c(-0.3, 0.3),
    expand = c(0, 0)) +
  coord_cartesian(clip = 'off')

ggsave(plot = p, paste0(PLT_DIR, '/Bias_9K_TEVp_ABC.png'),
  width = 3, height = 3)

p <- ggplot() +
  geom_logo(data = data_plot_lookup_MLD - 0.05, method = 'custom', seq_type = 'aa') +
  theme_SH() +
  scale_x_continuous('randomised position in TEVp',
    breaks = c(1, 2, 3), labels = 30:32) +
  scale_y_continuous('AA fraction',
    limits = c(-0.3, 0.3), 
    expand = c(0, 0)) +
  coord_cartesian(clip = 'off')

ggsave(plot = p, paste0(PLT_DIR, '/Bias_MLD_TEVp_ABC.png'),
  width = 3, height = 3)
