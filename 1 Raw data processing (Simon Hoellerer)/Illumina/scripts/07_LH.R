# import variables
ROOT_DIR <- Sys.getenv("ROOT_DIR")
OUT_DIR <- Sys.getenv("OUT_DIR")
RES_DIR <- Sys.getenv("RES_DIR")
PLT_DIR <- Sys.getenv("PLT_DIR")

# ROOT_DIR <- "/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/uASPIre_screen4_5/2022_LH_SH_S4_1"
# OUT_DIR <- "/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/uASPIre_screen4_5/2022_LH_SH_S4_1/output"
# RES_DIR <- "/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/uASPIre_screen4_5/2022_LH_SH_S4_1/results"

# change directory
setwd(ROOT_DIR)

# load libraries
library(tidyverse)
library(data.table)
library(matrixStats)

# define cutoff
CUTOFF <- 1
COL_NAMES <- c(paste0('N', 1:5), paste0('F', 1:5))

# read indices
IDX <- read.table(paste0(ROOT_DIR, '/Indices.txt'),
  header = T, sep = '\t')

# initialize dummy data
data <- data.frame(seq = 'X')

# read data into a list
for (L in c(1:5)) {
  for (R in c(3:4)) {
    for (state in c('N', 'F')) { # FF as F = FALSE
      # set input
      input_dat <- paste0(OUT_DIR, '/L', L, '_R', R, '_', state, '.txt')
      print(input_dat)
      
      # read data and index
      data_temp <- read.table(input_dat,
        colClasses = c('character', 'integer'), header = F, sep = '\t')
      names(data_temp) <- c('seq', paste0(state, L, R))

      # merge
      data <- data %>%
        full_join(., data_temp, by = 'seq')
    }
  }
}

# replace NA with 0
data[is.na(data)] <- 0

# calculate row sum
data$sum <- rowSums(data[, -1])
data <- data %>% arrange(desc(sum))

# split into libraries and samples
data_list <- list()

libs <- unique(IDX$lib)

for (idx_lib in libs) {
  # filter IDX list for each library
  IDX_lib <- IDX %>% filter(lib == idx_lib)

  # extract data
  data_temp <- data[, c('seq',
    paste0('N', IDX_lib$L, IDX_lib$R),
    paste0('F', IDX_lib$L, IDX_lib$R)
    )]
  names(data_temp) <- c('seq',
    paste0('N', IDX_lib$tp),
    paste0('F', IDX_lib$tp)
    )
  data_temp$sum <- rowSums(data_temp[, -1])

  # filter
  data_temp_HQ <- data_temp %>% filter(sum >= 1) %>% arrange(desc(sum))

  # add library
  data_temp_HQ$lib <- idx_lib

  # add to list
  data_list[[idx_lib]] <- data_temp_HQ
}

# combine
data_all <- as.data.frame(rbindlist(data_list))

# add full library
data_all$library <- NA
data_all$library[data_all$lib == 1] <- 'DMS'
data_all$library[data_all$lib == 2] <- 'MLD'

data <- data_all

# clean up and free some space
rm(data_all, data_list, data_temp_HQ, data_temp, L, R, state, IDX, IDX_lib, idx_lib,
  input_dat)
gc()

# find position of SpeI, HindIII and PstI
SpeI <- 'TCTAGT'
HindIII <- 'AAGCTT'
PstI <- 'CTGCAG'

pos_Spe <- as.data.frame(str_locate(data$seq, SpeI)) %>%
  rename(start_Spe = start, end_Spe = end)
pos_Hin <- as.data.frame(str_locate(data$seq, HindIII)) %>%
  rename(start_Hin = start, end_Hin = end)
pos_Pst <- as.data.frame(str_locate(data$seq, PstI)) %>%
  rename(start_Pst = start, end_Pst = end)

# combine
data <- cbind(data, pos_Spe, pos_Hin, pos_Pst)

# filter for correct cut sites
data$BC_TEVp <- substring(data$seq, data$end_Spe+1, data$start_Hin-1)
data$BC_TEVs <- substring(data$seq, data$end_Hin+1, data$start_Pst-1)

# crop
data <- data %>% select(BC_TEVp, BC_TEVs, all_of(COL_NAMES), library)

# split by libraries (for computational reasons)
DMS <- data %>% filter(library == 'DMS') %>% select(-library)
MLD <- data %>% filter(library == 'MLD') %>% select(-library)

# clean up and free some space
rm(data, pos_Spe, pos_Hin, pos_Pst, SpeI, HindIII, PstI)
gc()

# summarize by BCs
DMS <- DMS %>%
  group_by(BC_TEVp, BC_TEVs) %>%
  summarize_at(COL_NAMES, sum) %>%
  ungroup() %>%
  mutate(library = 'DMS')

MLD <- MLD %>%
  group_by(BC_TEVp, BC_TEVs) %>%
  summarize_at(COL_NAMES, sum) %>%
  ungroup() %>%
  mutate(library = 'MLD')

data <- rbind(DMS, MLD)

rm(DMS, MLD)
gc()

# recalculate row sum and row min
data$sum <- rowSums(data[, COL_NAMES])
data$min <- rowMins(as.matrix(data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))]))

# sort and remove missing barcodes
data <- data %>% arrange(desc(sum)) %>%
  filter(!is.na(BC_TEVp)) %>%
  filter(!is.na(BC_TEVs))

# write to file
write.table(data,
  file = paste0(RES_DIR, '/data_all.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')
