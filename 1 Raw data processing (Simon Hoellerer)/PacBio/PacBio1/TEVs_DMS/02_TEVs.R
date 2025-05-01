# set variables
ROOT_DIR <- '/links/groups/panke/07_Alumni/PhD/40_2023_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio1/TEVs_DMS'
RAW_DIR <- paste0(ROOT_DIR, '/raw_data')
TMP_DIR <- paste0(ROOT_DIR, '/tmp')
OUT_DIR <- paste0(ROOT_DIR, '/output')
RES_DIR <- paste0(ROOT_DIR, '/results')
PLT_DIR <- paste0(ROOT_DIR, '/plots')

# change directory
setwd(OUT_DIR)

# load libraries
library(tidyverse)
library(ggseqlogo)
library(spgs)
library(bioseq)
library(stringdist)

# read data
data_1 <- read.table('TEVs_left.txt',
  header = F, colClasses = c(rep('integer', 3), 'character'), sep = '\t')
data_2 <- read.table('TEVs_right.txt',
  header = F, colClasses = c(rep('integer', 3), 'character'), sep = '\t')
data_3 <- read.table('BC_left.txt',
  header = F, colClasses = c(rep('integer', 3), 'character'), sep = '\t')
data_4 <- read.table('BC_right.txt',
  header = F, colClasses = c(rep('integer', 3), 'character'), sep = '\t')

names(data_1) <- c('line', 'TLS', 'TLE', 'seq')
names(data_2) <- c('line', 'TRS', 'TRE', 'seq')
names(data_3) <- c('line', 'BLS', 'BLE', 'seq')
names(data_4) <- c('line', 'BRS', 'BRE', 'seq')

data_1 <- data_1 %>% select(line, seq, TLE)
data_2 <- data_2 %>% select(line, TRS)
data_3 <- data_3 %>% select(line, BLE)
data_4 <- data_4 %>% select(line, BRS)

# merge
data <- data_1 %>%
  inner_join(data_2, by = 'line') %>%
  inner_join(data_3, by = 'line') %>%
  inner_join(data_4, by = 'line')

# extract BC
data$BC <- substring(data$seq, data$BLE + 1, data$BRS)
data$TEVs <- substring(data$seq, data$TLE + 1, data$TRS)

# check length of BC and TEVs
temp <- nchar(data$BC) %>%
  table() %>%
  as.data.frame(stringsAsFactors = F) %>%
  arrange(desc(Freq))
names(temp) <- c('length_BC', 'freq')
class(temp$length_BC) <- 'integer'
temp <- temp %>% filter(length_BC < 20)

# plot
p <- ggplot(data = temp, aes(x = length_BC, y = freq)) +
  geom_bar(stat = 'identity', fill = 'dodgerblue3', colour = 'black', width = 0.80) +
  scale_x_continuous('barcode length (nt)',
    limits = c(7, 18), breaks = seq(7, 18, 1), expand = c(0, 0)) +
  scale_y_continuous('frequency',
    limits = c(0, 2500), expand = c(0, 0)) +
  theme_SH() +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/BC_TEVs_length.png'), plot = p,
  width = 3, height = 3, units = c('in'), scale = 1)

# check length of BC and TEVs
temp <- nchar(data$TEVs) %>%
  table() %>%
  as.data.frame(stringsAsFactors = F) %>%
  arrange(desc(Freq))
names(temp) <- c('length_TEVs', 'freq')
class(temp$length_TEVs) <- 'integer'

# plot
p <- ggplot(data = temp, aes(x = length_TEVs, y = freq)) +
  geom_bar(stat = 'identity', fill = 'dodgerblue3', colour = 'black', width = 0.80) +
  scale_x_continuous('TEVs length (nt)',
    limits = c(15, 25), breaks = seq(15, 25, 1), expand = c(0, 0)) +
  scale_y_continuous('frequency',
    limits = c(0, 2500), expand = c(0, 0)) +
  theme_SH() +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/TEVs_length.png'), plot = p,
  width = 3, height = 3, units = c('in'), scale = 1)

# filter
df <- data %>%
  filter(nchar(TEVs) == 21, nchar(BC) <= 12, nchar(BC) >= 8) %>%
  rename(BC_TEVs = BC)

df$trans <- seq_translate(dna(df$TEVs))

# write to file for processing
df_temp <- df %>%
  select(BC_TEVs, TEVs, trans, seq) %>%
  rename(TEVs_AA = trans)

# write to file
write.table(df_temp,
  file = paste0(RES_DIR, '/BC_TEVs_PacBio_with_reads.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

# plot BC
temp <- df %>% filter(nchar(BC_TEVs) == 10)
p <- ggplot() +
  geom_logo(data = temp$BC_TEVs, method = 'probability', seq_type = 'dna') +
  theme_SH() +
  scale_x_continuous('randomized position in TEVs BC',
    breaks = 1:10) +
  scale_y_continuous('base fraction', limits = c(0, 1)) +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/BC_TEVs_seq.png'), plot = p, width = 5,
  height = 2, units = c('in'), scale = 1)

# plot TEVs
p <- ggplot() +
  geom_logo(data = df$trans, method = 'probability', seq_type = 'aa') +
  theme_SH() +
  scale_x_continuous('randomized position in TEVs',
    breaks = seq(1, 7, 1)) +
  scale_y_continuous('AA fraction', limits = c(0, 1),
    breaks = seq(0, 1, 0.1)) +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/TEVs_seq.png'), plot = p, width = 5,
  height = 2, units = c('in'), scale = 1)

# find unique combinations
w <- df %>% group_by(BC_TEVs) %>% summarize(reads_BC_TEVs = n()) %>%
  ungroup() %>% arrange(desc(reads_BC_TEVs))

x <- df %>% group_by(BC_TEVs, TEVs) %>% summarize(reads_TEVsxBC = n()) %>%
  ungroup() %>% arrange(desc(reads_TEVsxBC))

# find number of unique TEVp per BC
y <- x %>% group_by(BC_TEVs) %>%
  summarize(TEVs_per_BC = n()) %>% ungroup() %>%
  arrange(desc(TEVs_per_BC))

# plot number of variants per BC
p <- ggplot(data = y, aes(x = TEVs_per_BC)) +
  geom_histogram(binwidth = 1, fill = 'dodgerblue3', color = 'black') +
  theme_SH() +
  scale_x_continuous('Number of TEVs per BC',
    limits = c(0, 5)) +
  scale_y_continuous('frequency',
    limits = c(0, 900), breaks = seq(0, 900, 100), expand = c(0,0)) +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/BC_TEVs_hist.png'), plot = p, width = 3,
  height = 3, units = c('in'), scale = 1)

z <- x %>% group_by(TEVs) %>%
  summarize(BCs_per_TEVs = n()) %>% ungroup() %>%
  arrange(desc(BCs_per_TEVs))

df_comb <- x %>%
  inner_join(w, by = 'BC_TEVs') %>%
  inner_join(y, by = 'BC_TEVs') %>%
  inner_join(z, by = 'TEVs')

# translate
df_comb$TEVs_trans <- seq_translate(dna(df_comb$TEVs))

# calculate fraction
df_comb <- df_comb %>% mutate(fraction_BC_TEVs = reads_TEVsxBC / reads_BC_TEVs)

# write to file
write.table(df_comb,
  file = paste0(RES_DIR, '/BC_TEVs_PacBio.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

# nrow(df_comb) # 1043
# length(unique(df_comb$BC_TEVs)) # 936
# length(unique(df_comb$TEVs_AA)) # 164
# stringdist('ENLYFQS', df_comb$TEVs_AA, method = 'hamming') %>% table()

# plot histogram
p <- ggplot(data = df_comb, aes(x = fraction_BC_TEVs)) +
  geom_histogram(binwidth = 0.05, fill = 'dodgerblue3', color = 'black') +
  theme_SH() +
  scale_x_continuous('Ratio of TEVs BC (white/blacklist)') +
  scale_y_continuous('frequency',
    limits = c(0, 900), breaks = seq(0, 900, 100), expand = c(0,0)) +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/ratio_hist.png'), plot = p, width = 3,
  height = 3, units = c('in'), scale = 1)

# done !
