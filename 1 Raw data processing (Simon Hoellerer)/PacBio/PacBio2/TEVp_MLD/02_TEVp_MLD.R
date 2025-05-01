# set variables
ROOT_DIR <- '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio2/TEVp_MLD'
RAW_DIR <- paste0(ROOT_DIR, '/raw_data')
TMP_DIR <- paste0(ROOT_DIR, '/tmp')
OUT_DIR <- paste0(ROOT_DIR, '/output')
RES_DIR <- paste0(ROOT_DIR, '/results')
PLT_DIR <- paste0(ROOT_DIR, '/plots')

# change directory
setwd(OUT_DIR)

# load libraries
library(tidyverse)
library(stringdist)
library(ggseqlogo)
library(bioseq)
library(spgs)

# read data
data <- read.table(paste0(TMP_DIR, '/03_CCS_MLD.txt'), sep = '\t', header = F,
  colClasses = 'character')
names(data) <- 'seq'

# add line number
data <- data %>% mutate(line = 1:nrow(.))

# read indices
i1 <- read.table('BC_L.txt', sep = '\t', header = F, colClasses = 'integer')
i2 <- read.table('BC_R.txt', sep = '\t', header = F, colClasses = 'integer')
i3 <- read.table('ABC_L.txt', sep = '\t', header = F, colClasses = 'integer')
i4 <- read.table('ABC_R.txt', sep = '\t', header = F, colClasses = 'integer')

i5 <- read.table('DEF_L.txt', sep = '\t', header = F, colClasses = 'integer')
i6 <- read.table('DEF_R.txt', sep = '\t', header = F, colClasses = 'integer')
i7 <- read.table('XYZ_L.txt', sep = '\t', header = F, colClasses = 'integer')
i8 <- read.table('XYZ_R.txt', sep = '\t', header = F, colClasses = 'integer')

names(i1) <- c('line', 'SBL', 'EBL')
names(i2) <- c('line', 'SBR', 'EBR')
names(i3) <- c('line', 'SAL', 'EAL')
names(i4) <- c('line', 'SAR', 'EAR')
names(i5) <- c('line', 'SDL', 'EDL')
names(i6) <- c('line', 'SDR', 'EDR')
names(i7) <- c('line', 'SXL', 'EXL')
names(i8) <- c('line', 'SXR', 'EXR')

# combine all
IDX <- i1 %>%
  inner_join(., i2, by = 'line') %>%
  inner_join(., i3, by = 'line') %>%
  inner_join(., i4, by = 'line') %>%
  inner_join(., i5, by = 'line') %>%
  inner_join(., i6, by = 'line') %>%
  inner_join(., i7, by = 'line') %>%
  inner_join(., i8, by = 'line')

# combine with data
data <- data %>% inner_join(., IDX, by = 'line')

# extract sequences
data <- data %>% mutate(BC_TEVp = substring(seq, EBL + 1, SBR))
data <- data %>% mutate(ABC = substring(seq, EAL + 1, SAR))
data <- data %>% mutate(DEF = substring(seq, EDL + 1, SDR))
data <- data %>% mutate(XYZ = substring(seq, EXL + 1, SXR))

# clean
data <- data %>% select(BC_TEVp, ABC, DEF, XYZ)
rm(i1, i2, i3, i4, i5, i6, i7, i8, IDX)

# count occurence
data <- data %>%
  group_by(BC_TEVp, ABC, DEF, XYZ) %>%
  summarize(reads = n()) %>%
  ungroup() %>%
  arrange(desc(reads))

# make barcode reverse complement (to fit plasmid orientation)
data <- data %>%
  mutate(BC_TEVp = reverseComplement(BC_TEVp, case = "upper"))

# plotting
for (current_position in c('BC_TEVp', 'ABC', 'DEF', 'XYZ')) {
  temp <- data %>%
    group_by_at(current_position) %>%
    summarize(sum = sum(reads)) %>%
    ungroup() %>%
    arrange(desc(sum)) %>%
    mutate(rank = rank(-sum))

  # plot
  p <- ggplot(temp, aes(x = log10(rank), y = log10(sum))) +
    geom_point(size = 0.1, color = 'dodgerblue3') +
    theme_SH() +
    coord_cartesian(clip = 'off')

  ggsave(p, file = paste0(PLT_DIR, '/Rank_', current_position, '.png'),
    width = 3, height = 3)
}

# add wt
ref_ABC <- 'ACCAGTCTG'
ref_DEF <- 'GATGGTCAG'
ref_XYZ <- 'TTTATGGTG'

# add variants that differ by only 1 from each set back to original
data$dist <- stringdist(data$ABC, ref_ABC, method = 'lv')
data$ABC[data$dist == 1] <- ref_ABC

data$dist <- stringdist(data$DEF, ref_DEF, method = 'lv')
data$DEF[data$dist == 1] <- ref_DEF

data$dist <- stringdist(data$XYZ, ref_XYZ, method = 'lv')
data$XYZ[data$dist == 1] <- ref_XYZ

data$dist <- NULL

# calculate reads again
data <- data %>% group_by(BC_TEVp, ABC, DEF, XYZ) %>%
  summarize(reads_TEVpxBC = sum(reads)) %>%
  ungroup() %>%
  group_by(BC_TEVp) %>%
  mutate(reads_BC_TEVp = sum(reads_TEVpxBC)) %>%
  mutate(TEVp_per_BC = n()) %>%
  ungroup() %>%
  group_by(ABC, DEF, XYZ) %>%
  mutate(reads_TEVp = sum(reads_TEVpxBC)) %>%
  mutate(BCs_per_TEVp = n()) %>%
  ungroup() %>%
  mutate(fraction_BC_TEVp = reads_TEVpxBC / reads_BC_TEVp) %>%
  arrange(desc(reads_TEVpxBC))

# add wildtype
data <- data %>%
  mutate(ABC_wt = (ABC == ref_ABC)) %>%
  mutate(DEF_wt = (DEF == ref_DEF)) %>%
  mutate(XYZ_wt = (XYZ == ref_XYZ))

### plots ###
temp <- data %>%
  filter(nchar(BC_TEVp) <= 20) %>% filter(nchar(BC_TEVp) >= 1) %>%
  filter(nchar(ABC) <= 20) %>% filter(nchar(ABC) >= 1) %>%
  filter(nchar(DEF) <= 20) %>% filter(nchar(DEF) >= 1) %>%
  filter(nchar(XYZ) <= 20) %>% filter(nchar(XYZ) >= 1)

p <- ggplot(data = temp, aes(x = nchar(BC_TEVp))) +
  geom_histogram(binwidth = 1, fill = 'dodgerblue3', color = 'black') +
  theme_SH() +
  scale_x_continuous('BC TEVp length') +
  scale_y_continuous('frequency') +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/BC_TEVp_length.png'), plot = p, width = 3,
  height = 3, units = c('in'), scale = 1)

# plot BCs
temp <- data %>% filter(nchar(BC_TEVp) == 15)
p <- ggplot() +
  geom_logo(data = temp$BC_TEVp, method = 'probability', seq_type = 'dna') +
  theme_SH() +
  scale_x_continuous('randomized position in TEVp BC',
    breaks = 1:15) +
  scale_y_continuous('base fraction',
    limits = c(0, 1)) +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/BC_TEVp_seq.png'), plot = p, width = 5,
  height = 2, units = c('in'), scale = 1)

# plot number of variants per BC
p <- data %>% group_by(ABC, DEF, XYZ, TEVp_per_BC) %>% summarize(n = n()) %>%
  ggplot(., aes(x = TEVp_per_BC)) +
    geom_histogram(binwidth = 1, fill = 'dodgerblue3', color = 'black') +
    theme_SH() +
    scale_x_continuous('Number of TEVp per BC') +
    scale_y_continuous('frequency') +
    coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/BC_TEVp_hist.png'), plot = p, width = 3,
  height = 3, units = c('in'), scale = 1)

p <- data %>% group_by(ABC, DEF, XYZ, BCs_per_TEVp) %>% summarize(n = n()) %>%
  ggplot(., aes(x = BCs_per_TEVp)) +
    geom_histogram(binwidth = 1, fill = 'dodgerblue3', color = 'black') +
    theme_SH() +
    scale_x_continuous('Number of BC per TEVp',
      limits = c(0, 20)) +
    scale_y_continuous('frequency') +
    coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/TEVp_BC_hist.png'), plot = p, width = 3,
  height = 3, units = c('in'), scale = 1)

# plot histogram
p <- ggplot(data = data, aes(x = fraction_BC_TEVp)) +
  geom_histogram(binwidth = 0.05, fill = 'dodgerblue3', color = 'black') +
  theme_SH() +
  scale_x_continuous('Ratio of TEVp BC (white/blacklist)') +
  scale_y_continuous('frequency') +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/ratio_hist.png'), plot = p, width = 3,
  height = 3)

# add libraries
df_WT <- data %>%
  filter(ABC_wt, DEF_wt, XYZ_wt) %>%
  mutate(library = 'WT')

df_ABC <- data %>%
  filter(!ABC_wt, DEF_wt, XYZ_wt) %>%
  mutate(library = 'ABC') %>%
  mutate(list = 'BL') %>%
  mutate(reason = 'wrong library')

df_DEF <- data %>%
  filter(ABC_wt, !DEF_wt, XYZ_wt) %>%
  mutate(library = 'DEF')

df_XYZ <- data %>%
  filter(ABC_wt, DEF_wt, !XYZ_wt) %>%
  mutate(library = 'XYZ')

df_ABCDEF <- data %>%
  filter(!ABC_wt, !DEF_wt, XYZ_wt) %>%
  mutate(library = 'ABCDEF')

df_ABCXYZ <- data %>%
  filter(!ABC_wt, DEF_wt, !XYZ_wt) %>%
  mutate(library = 'ABCXYZ')

df_DEFXYZ <- data %>%
  filter(ABC_wt, !DEF_wt, !XYZ_wt) %>%
  mutate(library = 'DEFXYZ') %>%
  mutate(list = 'BL') %>%
  mutate(reason = 'wrong library')

df_ABCDEFXYZ <- data %>%
  filter(!ABC_wt, !DEF_wt, !XYZ_wt) %>%
  mutate(library = 'ABCDEFXYZ')

# combine again
df_libs <- rbind(
  df_WT,
  df_DEF,
  df_XYZ,
  df_ABCDEF,
  df_ABCXYZ,
  df_ABCDEFXYZ)

### make black and whitelist
# blacklist all variants with wrong library lengths
BL_nchar <- df_libs %>%
  filter(!(nchar(ABC) == 9 & nchar(DEF) == 9 & nchar(XYZ) == 9)) %>%
  mutate(list = 'BL') %>%
  mutate(reason = 'frameshift')

WL <- df_libs %>%
  filter((nchar(ABC) == 9 & nchar(DEF) == 9 & nchar(XYZ) == 9)) %>%
  mutate(list = 'WL') %>%
  mutate(reason = 'WL')

# blacklist all with unclear BC assignment
BL_fraction <- WL %>%
  filter(fraction_BC_TEVp < 0.8) %>%
  mutate(list = 'BL') %>%
  mutate(reason = 'BL_WL_fraction')

# whitelist all with clear BC assignment
WL <- WL %>%
  filter(fraction_BC_TEVp >= 0.8) %>%
  mutate(list = 'WL') %>%
  mutate(reason = 'WL')

# add wrong libraries
BL <- rbind(BL_nchar, BL_fraction, df_ABC, df_DEFXYZ)

# remove from BL if present in WL
BL <- BL %>%
  mutate(isinwhite = BC_TEVp %in% WL$BC_TEVp) %>%
  filter(!isinwhite) %>%
  select(-isinwhite)

ref_TEVp_MLD <- rbind(WL, BL)

# translate
ref_TEVp_MLD <- ref_TEVp_MLD %>%
  mutate(TEVp_ABC = seq_translate(dna(ABC))) %>%
  mutate(TEVp_DEF = seq_translate(dna(DEF))) %>%
  mutate(TEVp_XYZ = seq_translate(dna(XYZ)))

# write to file
write.table(ref_TEVp_MLD,
  file = paste0(RES_DIR, '/BC_TEVp_PacBio_MLD.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

write.table(ref_TEVp_MLD,
  file = '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio_combined/TEVp/results/BC_TEVp_PacBio_MLD.txt',
  quote = F, row.names = F, col.names = T, sep = '\t')

# split again into variants that appear multiple times
ref_TEVp_MLD <- ref_TEVp_MLD %>%
  select(BC_TEVp, library, TEVp_ABC, TEVp_DEF, TEVp_XYZ, list, reason) %>%
  group_by(BC_TEVp) %>%
  mutate(n = n()) %>%
  ungroup()

ref_TEVp_MLD_1 <- ref_TEVp_MLD %>%
  filter(n == 1) %>%
  select(-n)

ref_TEVp_MLD_2 <- ref_TEVp_MLD %>%
  filter(n >= 2) %>%
  group_by(BC_TEVp) %>%
  summarize() %>%
  ungroup() %>%
  mutate(
    library = 'multiple',
    TEVp_ABC = 'XXX',
    TEVp_DEF = 'XXX',
    TEVp_XYZ = 'XXX',
    list = 'BL',
    reason = 'multiple')

# combine again
df_final <- rbind(
  ref_TEVp_MLD_1,
  ref_TEVp_MLD_2
  )

# write to file
write.table(df_final,
  file = paste0(RES_DIR, '/BC_TEVp_PacBio_MLD_clean.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

write.table(df_final,
  file = '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio_combined/TEVp/results/BC_TEVp_PacBio_MLD_clean.txt',
  quote = F, row.names = F, col.names = T, sep = '\t')

paste0('Final TEVp rows: ', nrow(df_final))
paste0('Final TEVp unique entries: ', length(unique(df_final$BC_TEVp)))

# done !




# # filter to remove very wrong reads
# data <- data %>%
#   filter(nchar(BC_TEVp) <= 20) %>% filter(nchar(BC_TEVp) >= 1) %>%
#   filter(nchar(ABC) <= 20) %>% filter(nchar(ABC) >= 1) %>%
#   filter(nchar(DEF) <= 20) %>% filter(nchar(DEF) >= 1) %>%
#   filter(nchar(XYZ) <= 20) %>% filter(nchar(XYZ) >= 1)
