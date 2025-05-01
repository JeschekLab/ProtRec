# set variables
ROOT_DIR <- '/links/groups/panke/07_Alumni/PhD/40_2023_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio2/TEVs_DMS'
RES_DIR <- paste0(ROOT_DIR, '/results')

# change directory
setwd(RES_DIR)

# load libraries
library(tidyverse)
library(seqinr)

# read data
data <- read.table(
  file = 'BC_TEVs_PacBio_with_reads.txt',
  header = T,
  colClasses = 'character')

# define function
extract <- function(BC){
  # check if BC is in file
  if(BC %in% data$BC_TEVs) {

    # filter data for requested barcode
    temp <- data %>%
      filter(BC_TEVs == BC)

    # calculate TEVs sequence and occurence
    temp2 <- temp %>%
      group_by(TEVs_AA) %>%
      summarize(n = n()) %>%
      arrange(desc(n))

    TEVs_AA <- temp2 %>% pull(TEVs_AA) %>% head(., 1)
    occurence <- temp2 %>% pull(n) %>% head(., 1)

    # print
    print(paste0('BC: ', BC))
    print(paste0('TEVs: ', TEVs_AA))
    print(paste0('N.: ', occurence))

    # write to file as fasta
    write.fasta(
      sequences = temp %>% pull(seq) %>% as.list(),
      names = paste0(BC, '_', 1:occurence),
      file.out = paste0(BC, '.fasta'))

    # print
    print(paste0('Written to file \'', BC, '.fasta\'!'))
  } else {
    print(paste0('ERROR: BC \'', BC, '\' not in file!'))
  }
}

# execute
extract('CGGTGTGGGG')

# done !
