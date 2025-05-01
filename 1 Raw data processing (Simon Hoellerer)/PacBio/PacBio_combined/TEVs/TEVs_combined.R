# set variables
ROOT_DIR <- '/links/groups/panke/07_Alumni/PhD/40_2023_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio_combined/TEVs'
PLT_DIR <- paste0(ROOT_DIR, '/plots')
RES_DIR <- paste0(ROOT_DIR, "/results")
PB_DIR1<- '/links/groups/panke/07_Alumni/PhD/40_2023_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio1'
PB_DIR2<- '/links/groups/panke/07_Alumni/PhD/40_2023_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio2'

# change directory
setwd(ROOT_DIR)

# load libraries
library(tidyverse)
library(VennDiagram)
library(ggseqlogo)
library(stringdist)

# create folders
dir.create(PLT_DIR)
dir.create(RES_DIR)

# read data
TEVs1 <- read.table(paste0(PB_DIR1, '/TEVs_DMS/results/BC_TEVs_PacBio.txt'),
  header = T, sep = '\t',
  colClasses = c(rep('character', 2),
    rep('integer', 4), 'character', 'numeric')) %>%
  mutate(run = 1) %>%
  rename(TEVs_DNA = TEVs, TEVs = TEVs_trans)

# read data
TEVs2 <- read.table(paste0(PB_DIR2, '/TEVs_DMS/results/BC_TEVs_PacBio.txt'),
  header = T, sep = '\t',
  colClasses = c(rep('character', 2),
    rep('integer', 4), 'character', 'numeric')) %>%
  mutate(run = 2) %>%
  rename(TEVs_DNA = TEVs, TEVs = TEVs_trans)

# Venn diagram
P1 <- TEVs1 %>% filter(fraction_BC_TEVs >= 0.8) %>% pull(BC_TEVs) %>% unique()
P2 <- TEVs2 %>% filter(fraction_BC_TEVs >= 0.8) %>% pull(BC_TEVs) %>% unique()

# make overlapping VennDiagram
venn.diagram(
  x = list(P1, P2),
  category.names = c('PacBio1' , 'PacBio2'),
  filename = paste0(PLT_DIR, '/BC_TEVs_Venn.pdf'),
  disable.logging = TRUE,
  OUTPUT = FALSE,
  fill = c('#2D4262', '#D09683'),
  alpha = c(0.5, 0.5),
  lwd = 0
)

# combine
data <- rbind(TEVs1, TEVs2)

# combine reads of identical barcodes
data <- data %>% group_by(BC_TEVs, TEVs, TEVs_DNA) %>%
  summarize(reads_TEVsxBC = sum(reads_TEVsxBC)) %>%
  ungroup() %>%
  group_by(BC_TEVs) %>%
  mutate(reads_BC_TEVs = sum(reads_TEVsxBC)) %>%
  mutate(TEVs_per_BC = n()) %>%
  ungroup() %>%
  group_by(TEVs) %>%
  mutate(reads_TEVs = sum(reads_TEVsxBC)) %>%
  mutate(BCs_per_TEVs = n()) %>%
  ungroup() %>%
  mutate(fraction_BC_TEVs = reads_TEVsxBC / reads_BC_TEVs) %>%
  arrange(desc(reads_TEVsxBC))

### plotting ###
p <- data %>% group_by(BC_TEVs) %>%
  summarize(length_BC = nchar(BC_TEVs)) %>%
  ggplot(., aes(x = length_BC)) +
    geom_histogram(binwidth = 1, fill = 'dodgerblue3', color = 'black') +
    scale_x_continuous('TEVs barcode length (nt)') +
    scale_y_continuous('frequency') +
    theme_SH() +
    coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/BC_TEVs_length.pdf'), plot = p,
  width = 3, height = 3)

# plot BC logo
temp <- data %>%
  group_by(BC_TEVs) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  mutate(length_BC = nchar(BC_TEVs)) %>%
  filter(length_BC == 10)

p <- ggplot() +
  geom_logo(data = temp$BC_TEVs, method = 'probability', seq_type = 'dna') +
  theme_SH() +
  scale_x_continuous('randomized position in TEVs BC',
    breaks = 1:10) +
  scale_y_continuous('base fraction', limits = c(0, 1)) +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/BC_TEVs_seq.pdf'), plot = p, width = 5,
  height = 2)

# plot TEVs as DNA
p <- ggplot() +
  geom_logo(data = data$TEVs_DNA, method = 'probability', seq_type = 'dna') +
  theme_SH() +
  scale_x_continuous('randomized position in TEVs',
    breaks = 1:21) +
  scale_y_continuous('base fraction', limits = c(0, 1)) +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/TEVs_DNA_seq.pdf'), plot = p, width = 5,
  height = 2)

# plot TEVs as AA
p <- ggplot() +
  geom_logo(data = data$TEVs, method = 'probability', seq_type = 'aa') +
  theme_SH() +
  scale_x_continuous('randomized position in TEVs',
    breaks = 1:7) +
  scale_y_continuous('AA fraction', limits = c(0, 1)) +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/TEVs_seq.pdf'), plot = p, width = 5,
  height = 2)

# find unique combinations
w <- data %>% group_by(BC_TEVs) %>% summarize(reads_BC_TEVs = n()) %>%
  ungroup() %>% arrange(desc(reads_BC_TEVs))

x <- data %>% group_by(BC_TEVs, TEVs) %>% summarize(reads_TEVsxBC = n()) %>%
  ungroup() %>% arrange(desc(reads_TEVsxBC))

# find number of unique TEVs per BC
y <- x %>% group_by(BC_TEVs) %>%
  summarize(TEVs_per_BC = n()) %>% ungroup() %>%
  arrange(desc(TEVs_per_BC))

# plot number of variants per BC
p <- ggplot(data = y %>% filter(TEVs_per_BC <= 10), aes(x = TEVs_per_BC)) +
  geom_histogram(binwidth = 1, fill = 'dodgerblue3', color = 'black') +
  theme_SH() +
  scale_x_continuous('Number of TEVs per BC') +
  scale_y_continuous('frequency') +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/TEVs_per_BC.pdf'), plot = p, width = 3,
  height = 3)

# add mutation
data <- data %>%
  mutate(n_mut = stringdist('ENLYFQS', TEVs, method = 'hamming'))

nrow(data) # 3110
length(unique(data$BC_TEVs)) # 2425
length(unique(data$TEVs_DNA)) # 257
table(data$n_mut)

# plot histogram
p <- ggplot(data = data, aes(x = fraction_BC_TEVs)) +
  geom_histogram(binwidth = 0.05, fill = 'dodgerblue3', color = 'black') +
  theme_SH() +
  scale_x_continuous('Ratio of TEVs BC (WL/BL)') +
  scale_y_continuous('frequency') +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/ratio_hist.pdf'), plot = p, width = 3,
  height = 3)

### make black and whitelist
# whitelist all with clear BC assignment
WL <- data %>%
  filter(fraction_BC_TEVs >= 0.8) %>%
  mutate(list = 'WL') %>%
  mutate(reason = 'WL')

# blacklist all with unclear BC assignment
BL <- data %>%
  filter(fraction_BC_TEVs < 0.8) %>%
  mutate(isinwhite = BC_TEVs %in% WL$BC_TEVs) %>%
  filter(!isinwhite) %>%
  mutate(list = 'BL') %>%
  mutate(reason = 'BL_WL_fraction') %>%
  select(-isinwhite)

df_combined <- rbind(WL, BL)

# write to file
write.table(df_combined,
  file = paste0(RES_DIR, '/BC_TEVs_PacBio_DMS.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

# split again into variants that appear multiple times
df_combined <- df_combined %>%
  group_by(BC_TEVs) %>%
  mutate(n = n()) %>%
  ungroup()

df_combined_1 <- df_combined %>%
  filter(n == 1) %>%
  select(BC_TEVs, TEVs, n_mut, list, reason)

df_combined_2 <- df_combined %>%
  filter(n >= 2) %>%
  group_by(BC_TEVs, list) %>%
  summarize() %>%
  ungroup() %>%
  mutate(
    TEVs = 'multiple',
    n_mut = 3,
    list = 'BL',
    reason = 'BL_WL_fraction')

df_final <- rbind(
  df_combined_1,
  df_combined_2
  )

write.table(df_final ,
  file = paste0(RES_DIR, '/BC_TEVs_PacBio_DMS_clean.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

paste0('Final TEVs rows: ', nrow(df_final))
paste0('Final TEVs unique entries: ', length(unique(df_final$BC_TEVs)))

# done !
