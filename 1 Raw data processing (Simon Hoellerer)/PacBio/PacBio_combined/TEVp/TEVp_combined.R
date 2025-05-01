# set variables
ROOT_DIR <- '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio_combined/TEVp'
PLT_DIR <- paste0(ROOT_DIR, '/plots')
RES_DIR <- paste0(ROOT_DIR, "/results")
PB_DIR1<- '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio1'
PB_DIR2<- '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio2'

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
TEVp1 <- read.table(paste0(PB_DIR1, '/TEVp_DMS/results/BC_TEVp_PacBio.txt'),
  header = T, sep = '\t') %>%
  mutate(run = 1)

TEVp2 <- read.table(paste0(PB_DIR2, '/TEVp_DMS/results/BC_TEVp_PacBio.txt'),
  header = T, sep = '\t') %>%
  mutate(run = 2)

# Venn diagram
P1 <- TEVp1 %>% filter(fraction_BC_TEVp >= 0.8) %>% pull(BC_TEVp) %>% unique()
P2 <- TEVp2 %>% filter(fraction_BC_TEVp >= 0.8) %>% pull(BC_TEVp) %>% unique()

# make overlapping VennDiagram
venn.diagram(
  x = list(P1, P2),
  category.names = c('PacBio1' , 'PacBio2'),
  filename = paste0(PLT_DIR, '/BC_TEVp_Venn.pdf'),
  disable.logging = TRUE,
  OUTPUT = FALSE,
  fill = c('#2D4262', '#D09683'),
  alpha = c(0.5, 0.5),
  lwd = 0
)

# combine
data <- rbind(TEVp1, TEVp2)

# count number of indels
indelstats <- data %>% group_by(indel, run) %>%
  summarize(reads = sum(reads_TEVpxBC)) %>%
  ungroup() %>%
  group_by(run) %>%
  mutate(total_reads = sum(reads)) %>%
  ungroup() %>%
  mutate(fraction = reads / total_reads)

# rename
indelstats$label <- 'correct'
indelstats$label[indelstats$indel] <- 'indel'

# plot
p <- ggplot(indelstats, aes(x = as.character(run), y = fraction, fill = label)) +
  geom_bar(position='stack', stat='identity') +
  theme_SH() +
  scale_x_discrete('PacBio run') +
  scale_y_continuous('fraction of reads', limits = c(0, 1),
    breaks = seq(0, 1, 0.1)) +
  scale_fill_manual(values = c('#2D4262', '#D09683')) +
  coord_cartesian(clip = 'off')

ggsave(plot = p, paste0(PLT_DIR, '/Indel_fraction.pdf'),
  width = 3, height = 3)

# combine reads of identical barcodes
data <- data %>% group_by(BC_TEVp, TEVp, n_mut, indel) %>%
  summarize(reads_TEVpxBC = sum(reads_TEVpxBC)) %>%
  ungroup() %>%
  group_by(BC_TEVp) %>%
  mutate(reads_BC_TEVp = sum(reads_TEVpxBC)) %>%
  mutate(TEVp_per_BC = n()) %>%
  ungroup() %>%
  group_by(TEVp) %>%
  mutate(reads_TEVp = sum(reads_TEVpxBC)) %>%
  mutate(BCs_per_TEVp = n()) %>%
  ungroup() %>%
  mutate(fraction_BC_TEVp = reads_TEVpxBC / reads_BC_TEVp) %>%
  arrange(desc(reads_TEVpxBC)) %>%
  mutate(comb = paste0(BC_TEVp, '_', TEVp))

### add variants back to non-indel if present in both
nrow(data) # 142249 reads
length(unique(data$TEVp)) # 39719 unique TEVp
length(unique(data$BC_TEVp)) # 105045 unique BCs

# split into with indel and without indel
df_noind <- data %>% filter(!indel)
df_indel <- data %>% filter(indel)

# add variants from indel to non-indel data set if they are present in both
df_indel$alsonoindel <- df_indel$comb %in% df_noind$comb
df_noind <- rbind(
  df_noind,
  df_indel %>% filter(alsonoindel) %>% mutate(indel = FALSE) %>% select(-alsonoindel)
  )

df_indel <- df_indel %>% filter(!alsonoindel) %>% select(-alsonoindel)

# now combine again
data <- rbind(df_noind, df_indel)

# calculate reads again
data <- data %>% group_by(BC_TEVp, TEVp, n_mut, indel) %>%
  summarize(reads_TEVpxBC = sum(reads_TEVpxBC)) %>%
  ungroup() %>%
  group_by(BC_TEVp) %>%
  mutate(reads_BC_TEVp = sum(reads_TEVpxBC)) %>%
  mutate(TEVp_per_BC = n()) %>%
  ungroup() %>%
  group_by(TEVp) %>%
  mutate(reads_TEVp = sum(reads_TEVpxBC)) %>%
  mutate(BCs_per_TEVp = n()) %>%
  ungroup() %>%
  mutate(fraction_BC_TEVp = reads_TEVpxBC / reads_BC_TEVp) %>%
  arrange(desc(reads_TEVpxBC))

### plotting ###
p <- data %>% group_by(BC_TEVp) %>%
  summarize(length_BC = nchar(BC_TEVp)) %>%
  ggplot(., aes(x = length_BC)) +
    geom_histogram(binwidth = 1, fill = 'dodgerblue3', color = 'black') +
    scale_x_continuous('TEVp barcode length (nt)') +
    scale_y_continuous('frequency') +
    theme_SH() +
    coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/BC_TEVp_length.pdf'), plot = p,
  width = 3, height = 3)

# plot BC logo
temp <- data %>%
  group_by(BC_TEVp) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  mutate(length_BC = nchar(BC_TEVp)) %>%
  filter(length_BC == 15)

p <- ggplot() +
  geom_logo(data = temp$BC_TEVp, method = 'probability', seq_type = 'dna') +
  theme_SH() +
  scale_x_continuous('randomized position in TEVp BC',
    breaks = 1:15) +
  scale_y_continuous('base fraction', limits = c(0, 1)) +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/BC_TEVp_seq.pdf'), plot = p, width = 5,
  height = 2)

# find unique combinations
w <- data %>% group_by(BC_TEVp) %>% summarize(reads_BC_TEVp = n()) %>%
  ungroup() %>% arrange(desc(reads_BC_TEVp))

x <- data %>% group_by(BC_TEVp, TEVp) %>% summarize(reads_TEVpxBC = n()) %>%
  ungroup() %>% arrange(desc(reads_TEVpxBC))

# find number of unique TEVp per BC
y <- x %>% group_by(BC_TEVp) %>%
  summarize(TEVp_per_BC = n()) %>% ungroup() %>%
  arrange(desc(TEVp_per_BC))

# plot number of variants per BC
p <- ggplot(data = y %>% filter(TEVp_per_BC <= 10), aes(x = TEVp_per_BC)) +
  geom_histogram(binwidth = 1, fill = 'dodgerblue3', color = 'black') +
  theme_SH() +
  scale_x_continuous('Number of TEVp per BC') +
  scale_y_continuous('frequency') +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/TEVp_per_BC.pdf'), plot = p, width = 3,
  height = 3)

nrow(data) # 134528
length(unique(data$BC_TEVp)) # 105045
length(unique(data$TEVp)) # 39719
table(data$n_mut)

# plot histogram
p <- ggplot(data = data, aes(x = fraction_BC_TEVp)) +
  geom_histogram(binwidth = 0.05, fill = 'dodgerblue3', color = 'black') +
  theme_SH() +
  scale_x_continuous('Ratio of TEVp BC (WL/BL)') +
  scale_y_continuous('frequency') +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/ratio_hist.pdf'), plot = p, width = 3,
  height = 3)

# count
p <- data %>% group_by(n_mut) %>% summarise(n = sum(reads_TEVpxBC)) %>%
  mutate(frac = n / sum(n) * 100) %>%
  ggplot(data = ., aes(x = n_mut, y = frac)) +
    geom_bar(stat = 'identity', fill = 'dodgerblue3', colour = 'black') +
    scale_x_continuous('number of TEVp mutations') +
    scale_y_continuous('fraction (%)',
      limits = c(0, 100), breaks = seq(0, 100, 10)) +
    theme_SH() +
    coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/TEVp_mutations.pdf'), plot = p, width = 3,
  height = 3, units = c('in'), scale = 1)

### make black and whitelist
# whitelist all with clear BC assignment
WL <- data %>%
  filter(fraction_BC_TEVp >= 0.8) %>%
  mutate(list = 'WL') %>%
  mutate(reason = 'WL')

# blacklist all with unclear BC assignment
BL <- data %>%
  filter(fraction_BC_TEVp < 0.8) %>%
  mutate(isinwhite = BC_TEVp %in% WL$BC_TEVp) %>%
  filter(!isinwhite) %>%
  mutate(list = 'BL') %>%
  mutate(reason = 'BL_WL_fraction') %>%
  select(-isinwhite)

ref_TEVp_DMS <- rbind(WL, BL)

# correct DMS lookup table to remove X123X variants that are actually wildtype
# first split into single and multiple mutants
ref_0 <- ref_TEVp_DMS %>%
  filter(n_mut == 0)

ref_1 <- ref_TEVp_DMS %>%
  filter(n_mut == 1)

ref_2 <- ref_TEVp_DMS %>%
  filter(n_mut >= 2)

# remove synonymous mutants from single mutant data set
ref_1 <- ref_1 %>%
  mutate(AAs = gsub("[^a-zA-Z]", "", TEVp)) %>%
  mutate(original_AA = substring(AAs, 1, 1)) %>%
  mutate(mutated_AA = substring(AAs, 2, 2)) %>%
  mutate(same = (mutated_AA == original_AA))

# filter and clean-up
ref_1_WL <- ref_1 %>% filter(!same) %>%
  select(-AAs, -original_AA, -mutated_AA, -same)

ref_1_BL <- ref_1 %>% filter(same) %>%
  select(-AAs, -original_AA, -mutated_AA, -same) %>%
  mutate(list = 'BL') %>%
  mutate(reason = 'synonymous')

# combine again
df_combined <- rbind(
  ref_0,
  ref_1_WL,
  ref_1_BL,
  ref_2)

# write to file
write.table(df_combined,
  file = paste0(RES_DIR, '/BC_TEVp_PacBio_DMS.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

# split again into variants that appear multiple times
df_combined <- df_combined %>%
  group_by(BC_TEVp) %>%
  mutate(n = n()) %>%
  ungroup()

df_combined_1 <- df_combined %>%
  filter(n == 1) %>%
  select(BC_TEVp, TEVp, indel, n_mut, list, reason)

df_combined_2 <- df_combined %>%
  filter(n >= 2) %>%
  group_by(BC_TEVp, list) %>%
  summarize() %>%
  ungroup() %>%
  mutate(
    TEVp = 'multiple',
    indel = TRUE,
    n_mut = 3,
    list = 'BL',
    reason = 'BL_WL_fraction')

# combine again
df_final <- rbind(
  df_combined_1,
  df_combined_2
  )

write.table(df_final ,
  file = paste0(RES_DIR, '/BC_TEVp_PacBio_DMS_clean.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

paste0('Final TEVp rows: ', nrow(df_final))
paste0('Final TEVp unique entries: ', length(unique(df_final$BC_TEVp)))

# done !
