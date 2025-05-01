# import variables
ROOT_DIR <- "/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/uASPIre_screen4_5"

# change directory
setwd(ROOT_DIR)

# load libraries
library(tidyverse)

# define cutoff
CUTOFF <- 1
COL_NAMES <- c(paste0('N', 1:5), paste0('F', 1:5))

# initialize data
DMS <- data.frame()
MLD <- data.frame()

# read data
for (current_run in 1:3) {
  # read table
  temp <- read.table(
    file = paste0(ROOT_DIR, '/2022_LH_SH_S4_', current_run, '/results/data_all.txt'),
    header = T, sep = '\t',
    colClasses = c(rep('character', 2), rep('integer', 10),
      'character', rep('integer', 2)))
  # split into libraries
  temp_DMS <- temp %>% filter(library == 'DMS') %>% select(-library)
  temp_MLD <- temp %>% filter(library == 'MLD') %>% select(-library)
  # combine
  DMS <- rbind(DMS, temp_DMS)
  MLD <- rbind(MLD, temp_MLD)
  # clean up
  rm(temp, temp_DMS, temp_MLD)
}

# free some space
gc()

# summarize by BCs
DMS <- DMS %>% select(-sum, -min, -run)
MLD <- MLD %>% select(-sum, -min, -run)

DMS <- DMS %>%
  group_by(BC_TEVp, BC_TEVs) %>%
  summarize_at(COL_NAMES, sum) %>%
  ungroup() %>%
  mutate(library = 'DMS')

# write to file
write.table(DMS,
  file = paste0('data_DMS.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')
rm(DMS)
gc()

MLD <- MLD %>%
  group_by(BC_TEVp, BC_TEVs) %>%
  summarize_at(COL_NAMES, sum) %>%
  ungroup() %>%
  mutate(library = 'MLD')

# write to file
write.table(MLD,
  file = paste0('./2022_LH_SH_S4_combined/results/data_MLD.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')
rm(MLD)
gc()
