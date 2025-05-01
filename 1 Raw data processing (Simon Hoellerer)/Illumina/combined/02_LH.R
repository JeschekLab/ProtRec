# import variables
ROOT_DIR <- "/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/uASPIre_screen4_5/2022_LH_SH_S4_combined"
PLT_DIR <- paste0(ROOT_DIR, '/plots')
RES_DIR <- paste0(ROOT_DIR, '/results')

# change directory
setwd(ROOT_DIR)

# load libraries
library(tidyverse)
library(matrixStats)

# library(data.table)
# library(stringdist)
# library(ggseqlogo)
# library(scico)
# library(pals)
# library(stringr)
# library(spgs)
# library(VennDiagram)

# define cutoff
CUTOFF <- 1
COL_NAMES <- c(paste0('N', 1:5), paste0('F', 1:5))

# initialize data
data <- rbind(
  read.table(
    file = paste0(RES_DIR, '/data_DMS.txt'),
    header = T, sep = '\t',
    colClasses = c(rep('character', 2), rep('integer', 10), 'character')),
  read.table(
    file = paste0(RES_DIR, '/data_MLD.txt'),
    header = T, sep = '\t',
    colClasses = c(rep('character', 2), rep('integer', 10), 'character')))

# add row sum and row min
data$sum <- rowSums(data[, COL_NAMES])
data$min <- rowMins(as.matrix(data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))]))

# dot plot of TEVs BC over reads
tempS <- data %>% group_by(library, BC_TEVs) %>%
  summarize(total_reads = sum(sum)) %>%
  ungroup() %>%
  arrange(desc(total_reads)) %>%
  mutate(rank = 1:nrow(.))

# plot
setwd(PLT_DIR)

p <- ggplot(tempS, aes(x = rank, y = total_reads, color = library)) +
  geom_point(size = 0.1) +
  theme_SH() +
  scale_y_log10()

ggsave(p, file = paste0(PLT_DIR, '/rank.png'))








# read table
data_1 <- read.table('/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/results/00_data.txt',
  header = T, sep = '\t')

# some plots
p <- data %>% filter(library == 'DMS') %>% group_by(BC_TEVs) %>%
  summarize(n = n()) %>% arrange(desc(n)) %>% mutate(count = 1:nrow(.)) %>%
  mutate(library = 'DMS') %>%
  ggplot(., aes(x = count, y = n, color = library)) +
    geom_point(size = 0.2) +
    theme_SH() +
    scale_x_continuous('sorted TEVs barcodes',
      limits = c(0, 2000)) +
    scale_y_continuous('total reads',
      limits = c(0, 2500))

ggsave(paste0(PLOT_DIR, '/01_NGS_TEVs_BC_DMS.png'),
  plot = p, width = 4, height = 3, units = c('in'), scale = 1)

p <- data %>% filter(library == 'MLD') %>% group_by(BC_TEVs) %>%
  summarize(n = n()) %>% arrange(desc(n)) %>% mutate(count = 1:nrow(.)) %>%
  mutate(library = 'MLD') %>%
  ggplot(., aes(x = count, y = n, color = library)) +
    geom_point(size = 0.2) +
    theme_SH() +
    scale_x_continuous('sorted TEVs barcodes',
      limits = c(0, 200)) +
    scale_y_continuous('total reads',
      limits = c(0, 100000))

ggsave(paste0(PLOT_DIR, '/01_NGS_TEVs_BC_MLD.png'),
  plot = p, width = 4, height = 3, units = c('in'), scale = 1)


# should we do merging with reading errors? -> let's not for now

# read in BC list from Excel
BC_TEVs_manual <- read.table('/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/BC_TEVs.txt',
  header = T, sep = '\t') %>%
  mutate(source_TEVs = 'manual')
BC_control <- read.table('/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/BC_controls.txt',
  header = T, sep = '\t') %>%
  mutate(source_TEVp = 'control')
BC_validation <- read.table('/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/BC_validation.txt',
  header = T, sep = '\t') %>%
  mutate(source_TEVp = 'validation')
BC_TEVs_PacBio <- read.table('/home/hoelsimo/02_Sequencing/NGS/2021_LH_TEV_PacBio/BC_TEVs_PacBio.txt',
  header = T, sep = '\t') %>%
  mutate(source_TEVs = 'PacBio')
BC_TEVp_PacBio <- read.table('/home/hoelsimo/02_Sequencing/NGS/2023_LH_TEV_PacBio/BC_TEVp_PacBio.txt',
  header = T, sep = '\t') %>%
  mutate(source_TEVp = 'PacBio')

# merge with manual TEVs
data <- data %>% full_join(., BC_TEVs_manual, by = 'BC_TEVs')
data <- data %>% full_join(., BC_control, by = 'BC_TEVp')
data <- data %>% full_join(., BC_validation, by = 'BC_TEVp')

# # combine source
# data$source_TEVs_manual[is.na(data$source_TEVs_manual)] <- ''
# data$source_TEVp_control[is.na(data$source_TEVp_control)] <- ''
# data$source_TEVp_validation[is.na(data$source_TEVp_validation)] <- ''

# data <- data %>% mutate(source = paste0(source_TEVs_manual, source_TEVp_control, source_TEVp_validation))

################ DMS INTERNAL CONTROLS ################
# extract library
df_cont <- data %>%
  filter(library == 'DMS') %>%
  full_join(., BC_TEVs_manual, by = 'BC_TEVs') %>%
  full_join(., BC_control, by = 'BC_TEVp') %>%
  rename(source_TEVs = source_TEVs_manual) %>%
  rename(source_TEVp = source_TEVp_control) %>%
  rename(libraryID = libraryID_control) %>%
  rename(PlasmidID = PlasmidID_control) %>%
  rename(TEVp = source_TEVp_control) %>%

# extract controls
  

  rename(source_TEVs = source_TEVs_manual) %>%
  rename(source_TEVp = source_TEVp_control) %>%
  select(-libraryID_validation, -PlasmidID_validation, -source_TEVp_validation,,
    -source_TEVp_validation, -source)

# calculate fraction flipped
for (tp in c(1:5)) {
  df_cont[, paste0('t', tp)] <- df_cont[, paste0('F', tp)] / 
    (df_cont[, paste0('F', tp)] + df_cont[, paste0('N', tp)])
}

# calculate AUC
tps <- c(0, 60, 120, 180, 240)

# calculate AUC
x <- as.character(tps) # get time points
df_y <- df_cont[, paste0('t', 1:5)] %>% as.data.frame(.) # get y values
x <- as.integer(x)

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_y), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_y[,i]) + (x[i+1]-x[i])*(df_y[,i+1]-df_y[,i])*0.5
  m_auc[, i] <- auc
}

df_cont$AUC <- round(rowSums(m_auc / 240), 4)

# write to file
write.table(df_cont,
  file = '/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/results/01_DMS_internal_controls.txt',
  quote = F, row.names = F, col.names = T, sep = '\t')


################ DMS CONTROLS ################
# extract controls
df_cont <- data %>%
  filter(library == 'DMS') %>%
  filter(source_TEVp_control == 'control') %>%
  rename(source_TEVp = source_TEVp_control) %>%
  rename(libraryID = libraryID_control) %>%
  rename(PlasmidID = PlasmidID_control) %>%
  rename(TEVp = source_TEVp_control) %>%
  select(-libraryID_validation, -PlasmidID_validation, -source_TEVp_validation,
    -source_TEVp_validation, -source, -TEVs_manual, -TEVs_trans)

# merge with TEVs from PacBio
BC_TEVs_PacBio_filtered <- BC_TEVs_PacBio %>% filter(fraction_BC_TEVs >= 0.9)
df_cont <- df_cont %>%
  full_join(., BC_TEVs_PacBio_filtered, by = 'BC_TEVs') %>%
  rename(source_TEVs = TEVs_PacBio)


# table(df_cont$TEVs_trans, useNA = 'always') %>% sort()

# merge reads with the same TEVs_trans
df_cont <- df_cont %>%
  group_by(TEVp, TEVs_trans, library, PlasmidID, libraryID, source_TEVp, Note,
    source_TEVs) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5)) %>%
  ungroup(); dim(df_cont)

df_cont$min <- rowMins(as.matrix(df_cont[, c(paste0('N', 1:5))] + df_cont[, c(paste0('F', 1:5))]))
df_cont$sum <- rowSums(as.matrix(df_cont[, c(paste0('N', 1:5), paste0('F', 1:5))]))

# calculate fraction flipped
for (tp in c(1:5)) {
  df_cont[, paste0('t', tp)] <- df_cont[, paste0('F', tp)] / 
    (df_cont[, paste0('F', tp)] + df_cont[, paste0('N', tp)])
}

# filter
data_HQ <- df_cont %>% filter(min >= 1) %>% arrange(desc(sum))

# calculate AUC
tps <- c(0, 60, 120, 180, 240)

# calculate AUC
x <- as.character(tps) # get time points
df_y <- data_HQ[, paste0('t', 1:5)] %>% as.data.frame(.) # get y values
x <- as.integer(x)

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_y), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_y[,i]) + (x[i+1]-x[i])*(df_y[,i+1]-df_y[,i])*0.5
  m_auc[, i] <- auc
}

data_HQ$AUC <- round(rowSums(m_auc/240), 4)

write.table(data_HQ, file = '01_DMS_controls.txt',
  quote = F, row.names = F, col.names = T, sep = '\t')


################ DMS VALIDATION MANUAL ################
# extract library
df_cont <- data %>%
  filter(library == 'DMS') %>%
  inner_join(., BC_TEVs_manual, by = 'BC_TEVs') %>%
  inner_join(., BC_validation, by = 'BC_TEVp') %>%
  rename(libraryID = libraryID_validation) %>%
  rename(PlasmidID = PlasmidID_validation)

# calculate fraction flipped
for (tp in c(1:5)) {
  df_cont[, paste0('t', tp)] <- df_cont[, paste0('F', tp)] / 
    (df_cont[, paste0('F', tp)] + df_cont[, paste0('N', tp)])
}

# calculate AUC
tps <- c(0, 60, 120, 180, 240)

# calculate AUC
x <- as.character(tps) # get time points
df_y <- df_cont[, paste0('t', 1:5)] %>% as.data.frame(.) # get y values
x <- as.integer(x)

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_y), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_y[,i]) + (x[i+1]-x[i])*(df_y[,i+1]-df_y[,i])*0.5
  m_auc[, i] <- auc
}

df_cont$AUC <- round(rowSums(m_auc / 240), 4)

# write to file
write.table(df_cont,
  file = '/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/results/01_DMS_validation_manual.txt',
  quote = F, row.names = F, col.names = T, sep = '\t')

################ DMS VALIDATION PacBio ################
# extract controls
df_cont <- data %>%
  filter(library == 'DMS') %>%
  filter(source_TEVp_validation == 'validation') %>%
  rename(libraryID = libraryID_validation) %>%
  rename(PlasmidID = PlasmidID_validation) %>%
  rename(TEVp = source_TEVp_validation) %>%
  rename(source_TEVp = source_TEVp_validation) %>%
  select(-libraryID_control, -PlasmidID_control, -source_TEVp_control,
    -source_TEVp_control, -source, -TEVs_manual, -TEVs_trans)

# merge with TEVs from PacBio
BC_TEVs_PacBio_filtered <- BC_TEVs_PacBio %>% filter(fraction_BC_TEVs >= 0.9)
df_cont <- df_cont %>%
  full_join(., BC_TEVs_PacBio_filtered, by = 'BC_TEVs') %>%
  rename(source_TEVs = TEVs_PacBio)


# table(df_cont$TEVs_trans, useNA = 'always') %>% sort()

# merge reads with the same TEVs_trans
df_cont <- df_cont %>%
  group_by(TEVp, TEVs_trans, library, PlasmidID, libraryID, source_TEVp, Note,
    source_TEVs) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5)) %>%
  ungroup(); dim(df_cont)

df_cont$min <- rowMins(as.matrix(df_cont[, c(paste0('N', 1:5))] + df_cont[, c(paste0('F', 1:5))]))
df_cont$sum <- rowSums(as.matrix(df_cont[, c(paste0('N', 1:5), paste0('F', 1:5))]))

# calculate fraction flipped
for (tp in c(1:5)) {
  df_cont[, paste0('t', tp)] <- df_cont[, paste0('F', tp)] / 
    (df_cont[, paste0('F', tp)] + df_cont[, paste0('N', tp)])
}

# filter
data_HQ <- df_cont %>% filter(min >= 1) %>% arrange(desc(sum))

# calculate AUC
tps <- c(0, 60, 120, 180, 240)

# calculate AUC
x <- as.character(tps) # get time points
df_y <- data_HQ[, paste0('t', 1:5)] %>% as.data.frame(.) # get y values
x <- as.integer(x)

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_y), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_y[,i]) + (x[i+1]-x[i])*(df_y[,i+1]-df_y[,i])*0.5
  m_auc[, i] <- auc
}

data_HQ$AUC <- round(rowSums(m_auc/240), 4)

write.table(data_HQ, file = '01_DMS_validation_PacBio.txt',
  quote = F, row.names = F, col.names = T, sep = '\t')

df_cont <- data %>%

  rename(source_TEVs = source_TEVs_manual) %>%
  rename(source_TEVp = source_TEVp_control) %>%
  rename(libraryID = libraryID_control) %>%
  rename(PlasmidID = PlasmidID_control) %>%
  rename(TEVp = source_TEVp_control) %>%


################ DMS DATA ################
# first filter PacBio data
BC_TEVs_PacBio_filtered <- BC_TEVs_PacBio %>% filter(fraction_BC_TEVs >= 0.9)
BC_TEVp_PacBio_filtered <- BC_TEVp_PacBio %>% filter(fraction_BC_TEVp >= 0.9)

data_DMS <- data %>%
  filter(library == 'DMS') %>%
  full_join(., BC_TEVs_PacBio_filtered, by = 'BC_TEVs') %>%
  full_join(., BC_TEVp_PacBio_filtered, by = 'BC_TEVp')

# merge reads with the same TEVs_trans and TEVp
data_sum <- data_DMS %>%
  group_by(TEVp, TEVs_trans, library, source_TEVp, source_TEVs, n_mut) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5)) %>%
  ungroup(); dim(data_sum)

data_sum$min <- rowMins(as.matrix(data_sum[, c(paste0('N', 1:5))] + data_sum[, c(paste0('F', 1:5))]))
data_sum$sum <- rowSums(as.matrix(data_sum[, c(paste0('N', 1:5), paste0('F', 1:5))]))

# calculate fraction flipped
for (tp in c(1:5)) {
  data_sum[, paste0('t', tp)] <- data_sum[, paste0('F', tp)] / 
    (data_sum[, paste0('F', tp)] + data_sum[, paste0('N', tp)])
}

# filter
data_HQ <- data_sum %>% filter(min >= 1) %>% arrange(desc(sum))

# calculate AUC
tps <- c(0, 60, 120, 180, 240)

# calculate AUC
x <- as.character(tps) # get time points
df_y <- data_HQ[, paste0('t', 1:5)] %>% as.data.frame(.) # get y values
x <- as.integer(x)

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_y), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_y[,i]) + (x[i+1]-x[i])*(df_y[,i+1]-df_y[,i])*0.5
  m_auc[, i] <- auc
}

data_HQ$AUC <- round(rowSums(m_auc/240), 4)

data_HQ %>% group_by(source_TEVs) %>% summarize(n = n())
data_HQ %>% group_by(source_TEVs) %>% summarize(sum = sum(sum))
data_HQ %>% group_by(source_TEVp) %>% summarize(n = n())
data_HQ %>% group_by(source_TEVp) %>% summarize(sum = sum(sum))


write.table(data_HQ, file = '01_DMS_validation_PacBio.txt',
  quote = F, row.names = F, col.names = T, sep = '\t')

################ MLD ################
df_MLD <- data %>%
  filter(library == 'MLD') # %>%
  # left_join(., BC_TEVp_PacBio, by = 'BC_TEVp') ### DATA MISSING!!!

# map back TEVs - first step find unique TEVs
MAX_ERROR <- 1

df_correct <- df_MLD %>% filter(BC_TEVs %in% BC_TEVs_manual$BC_TEVs)
df_unmapped <- df_MLD %>% filter(!(BC_TEVs %in% BC_TEVs_manual$BC_TEVs))

unmapped_BCs <- df_unmapped %>% group_by(BC_TEVs) %>%
  summarize(reads = sum(sum)) %>% arrange(desc(reads)) %>% ungroup()

M <- stringdistmatrix(unmapped_BCs$BC_TEVs, BC_TEVs_manual$BC_TEVs,
  method = 'lv') %>% as.data.frame(., stringsAsFactors = F)
names(M) <- BC_TEVs_manual$BC_TEVs
M$BC_TEVs <- unmapped_BCs$BC_TEVs

# find the minimum distance per row
M$min <- rowMins(as.matrix(M[, 1:21]))

# find which column contains the minimum
M$which_min <- max.col(-M[, 1:21])

# count the appearance of the minimum
M$n_min <- rowSums(M[, 1:21] == M$min)

# make printable column
M$best_match <- names(M)[M$which_min]

# ignore lines that have their minimum more than once
M$best_match[M$n_min != 1] <- NA

# ignore lines that have a minimum distance greater than 1
M$best_match[M$min > MAX_ERROR] <- NA

# extract mapping columns
M <- M %>% select(BC_TEVs, best_match)

# merge back
df_mapped <- df_unmapped %>% left_join(., M, by = 'BC_TEVs')

# map BCs
df_correct <- df_correct %>%
  mutate(original = BC_TEVs) %>%
  mutate(mapped = FALSE) %>%
  left_join(., BC_TEVs_manual, by = 'BC_TEVs')

df_unmapped <- df_mapped %>%
  filter(is.na(best_match)) %>%
  mutate(original = BC_TEVs) %>%
  mutate(mapped = FALSE) %>%
  select(-best_match) %>%
  left_join(., BC_TEVs_manual, by = 'BC_TEVs')

df_mapped <- df_mapped %>%
  filter(!is.na(best_match)) %>%
  mutate(original = BC_TEVs) %>%
  mutate(mapped = TRUE) %>%
  mutate(BC_TEVs = best_match) %>%
  select(-best_match) %>%
  left_join(., BC_TEVs_manual, by = 'BC_TEVs')

# combine data frames
df_MLD_mapped <- rbind(df_correct, df_mapped, df_unmapped)

p <- df_MLD_mapped %>% group_by(TEVs_trans) %>% summarize(n = n()) %>%
  arrange(desc(n)) %>% mutate(fraction = n / sum(n)) %>%
  mutate(P1prime = substring(TEVs_trans, 7, 7)) %>%
  ggplot(., aes(x = P1prime, y = fraction)) +
  geom_bar(position = 'dodge', stat = 'identity',
    color = 'black', fill = 'dodgerblue3') +
  theme_SH() +
  scale_y_continuous('fraction of all reads',
    limits = c(0, 0.25), expand = c(0, 0))

ggsave(paste0(PLOT_DIR, '/02_MLD_P1_mapped.png'),
  plot = p, width = 4, height = 3, units = c('in'), scale = 1)

p <- df_MLD_mapped %>% group_by(BC_TEVs) %>%
  summarize(n = n()) %>% arrange(desc(n)) %>% mutate(count = 1:nrow(.)) %>%
  mutate(library = 'MLD') %>%
  ggplot(., aes(x = count, y = n, color = library)) +
    geom_point(size = 0.2) +
    theme_SH() +
    scale_x_continuous('sorted TEVs barcodes',
      limits = c(0, 200)) +
    scale_y_continuous('total reads',
      limits = c(0, 125000), breaks = seq(0, 125000, 25000))

ggsave(paste0(PLOT_DIR, '/01_NGS_TEVs_BC_MLD_mapped.png'),
  plot = p, width = 4, height = 3, units = c('in'), scale = 1)


# temp <- df_MLD_mapped %>% group_by(BC_TEVp) %>% summarize(sum = sum(sum)) %>%
#   arrange(desc(sum)) %>% mutate(count = 1:nrow(.))

# p <- ggplot(temp, aes(x = count, y = sum)) +
#   geom_point(size = 0.2, color = 'dodgerblue3') +
#   theme_SH() +
#   scale_x_continuous('sorted TEVp barcodes',
#     limits = c(0, 200000)) +
#   scale_y_continuous('total reads',

# ggsave(paste0(PLOT_DIR, '/01_NGS_TEVp_BC_MLD.png'),
#   plot = p, width = 3, height = 3, units = c('in'), scale = 1)



































# how many BCs in Illumina and PacBio
a <- BC_TEVs_PacBio$BC_TEVs %>% unique()
c <- df_cont$BC_TEVs %>% unique()

overlap <- data.frame(
  PacBio = c(sum(a %in% a), sum(a %in% c)),
  Illumina = c(sum(c %in% a), sum(c %in% c)))
rownames(overlap) <- c('PacBio', 'Illumina')
overlap %>% kable()

# make overlapping VennDiagram
myCol <- brewer.pal(3, 'Pastel2')

venn.diagram(
  x = list(a, b, c),
  category.names = c('PacBio' , 'validation' , 'controls'),
  filename = paste0(PLOT_DIR, '/01_Venn.png'),
  output = FALSE,
  disable.logging = TRUE,

  # Output features
  imagetype='png' ,
  height = 960 , 
  width = 960 , 
  resolution = 300,
  compression = 'lzw',
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = 'bold',
  fontfamily = 'sans',
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = 'bold',
  cat.default.pos = 'text',
  cat.fontfamily = 'sans',
  rotation = 1
)





















################ MLD INTERNAL CONTROLS ################
# extract controls
df_cont <- data %>%
  filter(library == 'MLD') %>%
  filter(source_TEVp_control == 'control') %>%
  filter(TEVs_manual == 'manual') %>%
  rename(libraryID = libraryID_control) %>%
  rename(PlasmidID = PlasmidID_control) %>%
  rename(TEVp = source_TEVp_control) %>%
  rename(source_TEVs = TEVs_manual) %>%
  rename(source_TEVp = source_TEVp_control) %>%
  select(-libraryID_validation, -PlasmidID_validation, -source_TEVp_validation,
    -source_TEVp_validation, -source)

# calculate fraction flipped
for (tp in c(1:5)) {
  df_cont[, paste0('t', tp)] <- df_cont[, paste0('F', tp)] / 
    (df_cont[, paste0('F', tp)] + df_cont[, paste0('N', tp)])
}

# calculate AUC
tps <- c(0, 60, 120, 180, 240)

# calculate AUC
x <- as.character(tps) # get time points
df_y <- df_cont[, paste0('t', 1:5)] %>% as.data.frame(.) # get y values
x <- as.integer(x)

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_y), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_y[,i]) + (x[i+1]-x[i])*(df_y[,i+1]-df_y[,i])*0.5
  m_auc[, i] <- auc
}

df_cont$AUC <- round(rowSums(m_auc / 240), 4)

# write to file
write.table(df_cont,
  file = '/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/results/01_MLD_internal_controls.txt',
  quote = F, row.names = F, col.names = T, sep = '\t')




























































# # rename
# df_cont <- df_cont %>%
#   rename(source_TEVs = TEVs_PacBio) %>%
#   rename(TEVs_trans = TEVs_AA)


#    %>%
#   filter(TEVs_manual == 'manual') %>%
  
  
#   rename(source_TEVs = TEVs_manual) %>%
  
  

# # calculate fraction flipped
# for (tp in c(1:5)) {
#   df_cont[, paste0('t', tp)] <- df_cont[, paste0('F', tp)] / 
#     (df_cont[, paste0('F', tp)] + df_cont[, paste0('N', tp)])
# }

# # calculate AUC
# tps <- c(0, 60, 120, 180, 240)

# # calculate AUC
# x <- as.character(tps) # get time points
# df_y <- df_cont[, paste0('t', 1:5)] %>% as.data.frame(.) # get y values
# x <- as.integer(x)

# # initialize helper matrix for AUC calculation
# m_auc <- matrix(nrow = nrow(df_y), ncol = (length(x)-1))

# # loop through columns and calculate AUC
# for (i in 1:(length(x)-1)) {
#   print(i)
#   auc <- (x[i+1]-x[i])*(df_y[,i]) + (x[i+1]-x[i])*(df_y[,i+1]-df_y[,i])*0.5
#   m_auc[, i] <- auc
# }

# df_cont$AUC <- round(rowSums(m_auc / 240), 4)

# # write to file
# write.table(df_cont,
#   file = '/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/results/01_DMS_controls.txt',
#   quote = F, row.names = F, col.names = T, sep = '\t')















# merge




data_DMS <- data_DMS %>% group_by(mutation_TEVp, trans_TEVs) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5))

data_DMS$min <- rowMins(as.matrix(data_DMS[, c(paste0('N', 1:5))] + data_DMS[, c(paste0('F', 1:5))]))
data_DMS$sum <- rowSums(as.matrix(data_DMS[, c(paste0('N', 1:5), paste0('F', 1:5))]))

# calculate fraction flipped
for (tp in c(1:5)) {
  data_DMS[, paste0('t', tp)] <- data_DMS[, paste0('F', tp)] / 
    (data_DMS[, paste0('F', tp)] + data_DMS[, paste0('N', tp)])
}

# filter
data_HQ <- data_DMS %>% filter(min >= 10) %>% arrange(desc(sum))

# calculate AUC
tps <- c(0, 60, 120, 180, 240)

# calculate AUC
x <- as.character(tps) # get time points
df_y <- data_HQ[, paste0('t', 1:5)] %>% as.data.frame(.) # get y values
x <- as.integer(x)

# initialize helper matrix for AUC calculation
m_auc <- matrix(nrow = nrow(df_y), ncol = (length(x)-1))

# loop through columns and calculate AUC
for (i in 1:(length(x)-1)) {
  print(i)
  auc <- (x[i+1]-x[i])*(df_y[,i]) + (x[i+1]-x[i])*(df_y[,i+1]-df_y[,i])*0.5
  m_auc[, i] <- auc
}

data_HQ$AUC <- round(rowSums(m_auc), 4)

# line plot
n <- 1000
subsample <- data_HQ[1:n, c(paste0('t', 1:5), 'AUC')] %>% arrange(AUC)

# make matrix
M <- as.matrix(subsample[, 1:6])
rownames(M) <- 1:n

# adapt matrix for diverging color
n_points <- length(seq(0, 240, 5))
plot_data <- as.data.frame(matrix(nrow = n, ncol = n_points),
  stringsAsFactors = F)

# change col and rownames
colnames(plot_data) <- seq(0, 240, 5)
rownames(plot_data) <- 1:n

# vectors
y_known <- tps
x_known <- (y_known / 5) + 1

for (i in 1:4) {
  x1 <- x_known[i]
  x2 <- x_known[i + 1]
  data_x1 <- M[, i]
  data_x2 <- M[, i + 1]
  len <- x2 - x1 + 1
  for (j in 1:n) {
    plot_data[j,x1:x2] <- seq(data_x1[j], data_x2[j], length = len)
  }
}

plot_data$ID <- 1:nrow(plot_data)

plot_melted <- reshape2::melt(plot_data, id = 'ID')
names(plot_melted) <- c('ID', 'tp', 'fraction')

p_line <- ggplot(plot_melted, aes(x = tp, y = fraction, group = ID)) + 
  geom_line(col = 'black', alpha = 0.5)+
  scale_colour_scico(palette = 'batlow', limits = c(0, 1)) +
  scale_x_discrete('time after induction (h)',
    breaks = seq(0, 240, 60), labels = c(0, 1, 2, 3, 4), expand = c(0,0)) +
  scale_y_continuous('fraction flipped (%)',
    breaks = seq(0, 1, 0.25), limits = c(0, 1), expand = c(0,0)) +
  theme_SH() +
  theme(legend.position = 'none')

ggsave(paste0('01_lineplot_lib_', n, '.png'),
  plot = p_line, width = 5, height = 5, units = c('in'), scale = 1)

# 
data_HQ %>% select(mutation_TEVp, trans_TEVs, min, sum, AUC, t1) %>%
  rename(TEVp = mutation_TEVp, TEVs = trans_TEVs, t0 = t1) %>%
  mutate(AUC = round(AUC / 240, 4), t0 = round(t0, 3)) %>%
  write.table(., file = '01_data_DMS.txt',
    quote = F, row.names = F, col.names = T, sep = '\t')


# check some varinats
data_HQ %>% select(mutation_TEVp, trans_TEVs, min, sum, AUC, t1) %>%
  arrange(desc(AUC)) %>% as.data.frame() %>% head(., 30)







data_sum <- data %>%
  group_by(BC_TEVp, BC_TEVs) %>%
  summarize(
    N13 = sum(N13), F13 = sum(F13),
    N14 = sum(N14), F14 = sum(F14),
    N23 = sum(N23), F23 = sum(F23),
    N24 = sum(N24), F24 = sum(F24),
    N33 = sum(N33), F33 = sum(F33),
    N34 = sum(N34), F34 = sum(F34),
    N43 = sum(N43), F43 = sum(F43),
    N44 = sum(N44), F44 = sum(F44),
    N53 = sum(N53), F53 = sum(F53),
    N54 = sum(N54), F54 = sum(F54)
    ) %>%
  ungroup()












