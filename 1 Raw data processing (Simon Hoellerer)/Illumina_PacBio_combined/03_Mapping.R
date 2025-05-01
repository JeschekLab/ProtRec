# import variables
ROOT_DIR <- '/links/groups/panke/07_Alumni/PhD/40_2023_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/uASPIre_screen4_5/Illumina_PacBio_combined'
PLT_DIR <- paste0(ROOT_DIR, '/plots')
RES_DIR <- paste0(ROOT_DIR, '/results')
Illumina_DIR <- '/links/groups/panke/07_Alumni/PhD/40_2023_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/uASPIre_screen4_5/2022_LH_SH_S4_combined'
PacBio_DIR <- '/links/groups/panke/07_Alumni/PhD/40_2023_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio_combined'

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



##################### DMS CONTROLS #####################
# extract library
data <- raw_DMS %>%
  inner_join(., ref_TEVp_controls, by = 'BC_TEVp') %>%
  inner_join(., ref_TEVs_DMS %>% filter(list_TEVs == 'WL'), by = 'BC_TEVs')

# write to file
write.table(data,
  file = paste0(RES_DIR, '/Data_DMS_controls_unmapped.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

print("Done!")



##################### DMS VARIANTS #####################
# extract library
data <- raw_DMS %>%
  inner_join(., ref_TEVp_DMS %>% filter(list_TEVp == 'WL'), by = 'BC_TEVp') %>%
  inner_join(., ref_TEVs_DMS %>% filter(list_TEVs == 'WL'), by = 'BC_TEVs')

# write to file
write.table(data,
  file = paste0(RES_DIR, '/Data_DMS_unmapped.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

print("Done!")

# # make new references
# ref_TEVp <- ref_TEVp_DMS %>% filter(BC_TEVp %in% unique(data$BC_TEVp))
# ref_TEVs <- ref_TEVs_DMS %>% filter(BC_TEVs %in% unique(data$BC_TEVs))

# # find wrong or missing TEVs
# missing_TEVs <- raw_DMS %>%
#   filter(!(BC_TEVs %in% ref_TEVs_DMS$BC_TEVs)) %>%
#   filter(!(BC_TEVs %in% ref_TEVs_controls$BC_TEVs)) %>%
#   filter(!(BC_TEVs %in% ref_TEVs_validation$BC_TEVs))

# # add length
# ref_TEVs <- ref_TEVs %>%
#   mutate(length = nchar(BC_TEVs))
# missing_TEVs <- missing_TEVs %>%
#   mutate(length = nchar(BC_TEVs))

# # filter for strange barcodes
# missing_TEVs <- missing_TEVs %>%
#   filter(nchar(BC_TEVs) %in% c(min(nchar(ref_TEVs$BC_TEVs)):max(nchar(ref_TEVs$BC_TEVs)))) %>%
#   arrange(desc(sum))

# # calculate distance to ref_TEVs variants
# # loop through lenghts
# LENGTHS <- sort(unique(ref_TEVs$length))

# # initialize mapped df
# mapped_df <- data.frame()

# for (current_length in LENGTHS) {
#   # make part of ref and map
#   print(paste0('Filtering length ', current_length, ' ...'))
#   ref_part <- ref_TEVs %>%
#     filter(length == current_length)
#   data_part <- missing_TEVs %>%
#     filter(length == current_length)

#   # extract unique BC_TEVs
#   data_part_uniq <- data_part %>%
#     pull(BC_TEVs) %>%
#     unique()

#   # calculate distance matrix
#   print(paste0('Calculate distance matrix of length ', current_length, ' ...'))
#   M <- stringdistmatrix(data_part_uniq, ref_part$BC_TEVs,
#     method = 'hamming', nthread = 11)

#   # find the minimum distance per row
#   print(paste0('Finding minimum of length ', current_length, ' ...'))
#   current_min <- rowMins(M)

#   # find which column contains the minimum
#   print(paste0('Finding which minimum of length ', current_length, ' ...'))
#   which_min <- max.col(-M)

#   # count the appearance of the minimum
#   print(paste0('Counting minimum of length ', current_length, ' ...'))
#   n_min <- rowSums(M == current_min)

#   # make printable column
#   print_names <- ref_part$BC_TEVs[which_min]

#   # ignore lines that have their minimum more than once
#   print_names[n_min != 1] <- NA

#   # ignore lines that have a minimum distance greater than 1
#   print_names[current_min > 1] <- NA

#   # add to data part
#   new_df <- data.frame(
#     BC_TEVs = data_part_uniq,
#     BC_TEVs_new = print_names)

#   # add to mapped M
#   mapped_df <- rbind(mapped_df, new_df)
# }

# # merge with missing_df
# mapped_DMS <- raw_DMS %>%
#   left_join(., mapped_df, by = 'BC_TEVs')

# # overwrite old data with new data
# mapped_DMS <- rbind(
#   mapped_DMS %>%
#     filter(is.na(BC_TEVs_new)),
#   mapped_DMS %>%
#     filter(!is.na(BC_TEVs_new)) %>%
#     mutate(BC_TEVs = BC_TEVs_new)
#   ) %>%
# select(-BC_TEVs_new)

# # combine reads and merge with data
# mapped_DMS <- mapped_DMS %>%
#   group_by(BC_TEVp, BC_TEVs) %>%
#   summarize(
#     N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
#     F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
#     sum = sum(sum)) %>%
#   ungroup()

# # find wrong or missing TEVp
# missing_TEVp <- mapped_DMS %>%
#   filter(!(BC_TEVp %in% ref_TEVp_DMS$BC_TEVp)) %>%
#   filter(!(BC_TEVp %in% ref_TEVp_controls$BC_TEVp)) %>%
#   filter(!(BC_TEVp %in% ref_TEVp_validation$BC_TEVp))

# # add length
# ref_TEVp <- ref_TEVp %>%
#   mutate(length = nchar(BC_TEVp))

# missing_TEVp <- missing_TEVp %>%
#   mutate(length = nchar(BC_TEVp))

# # filter for strange barcodes
# missing_TEVp <- missing_TEVp %>%
#   filter(nchar(BC_TEVp) %in% c(min(nchar(ref_TEVp$BC_TEVp)):max(nchar(ref_TEVp$BC_TEVp)))) %>%
#   arrange(desc(sum))

# # calculate distance to ref_TEVp variants
# # loop through lenghts
# LENGTHS <- sort(unique(ref_TEVp$length))

# # initialize mapped df
# mapped_df <- data.frame()

# for (current_length in LENGTHS) {
#   # make part of ref and map
#   print(paste0('Filtering length ', current_length, ' ...'))
#   ref_part <- ref_TEVp %>%
#     filter(length == current_length)
#   data_part <- missing_TEVp %>%
#     filter(length == current_length)

#   # extract unique BC_TEVp
#   data_part_uniq <- data_part %>%
#     pull(BC_TEVp) %>%
#     unique()

#   # calculate distance matrix
#   print(paste0('Calculate distance matrix of length ', current_length, ' ...'))
#   M <- stringdistmatrix(data_part_uniq, ref_part$BC_TEVp,
#     method = 'hamming', nthread = 11)

#   # find the minimum distance per row
#   print(paste0('Finding minimum of length ', current_length, ' ...'))
#   current_min <- rowMins(M)

#   # find which column contains the minimum
#   print(paste0('Finding which minimum of length ', current_length, ' ...'))
#   which_min <- max.col(-M)

#   # count the appearance of the minimum
#   print(paste0('Counting minimum of length ', current_length, ' ...'))
#   n_min <- rowSums(M == current_min)

#   # make printable column
#   print_names <- ref_part$BC_TEVp[which_min]

#   # ignore lines that have their minimum more than once
#   print_names[n_min != 1] <- NA

#   # ignore lines that have a minimum distance greater than 1
#   print_names[current_min > 1] <- NA

#   # add to data part
#   new_df <- data.frame(
#     BC_TEVp = data_part_uniq,
#     BC_TEVp_new = print_names)

#   # add to mapped M
#   mapped_df <- rbind(mapped_df, new_df)
# }

# # merge with missing_df
# mapped_DMS <- raw_DMS %>%
#   left_join(., mapped_df, by = 'BC_TEVp')

# # overwrite old data with new data
# mapped_DMS <- rbind(
#   mapped_DMS %>%
#     filter(is.na(BC_TEVp_new)),
#   mapped_DMS %>%
#     filter(!is.na(BC_TEVp_new)) %>%
#     mutate(BC_TEVp = BC_TEVp_new)
#   ) %>%
# select(-BC_TEVp_new)

# # combine reads and merge with data
# mapped_DMS <- mapped_DMS %>%
#   group_by(BC_TEVp, BC_TEVp) %>%
#   summarize(
#     N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
#     F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
#     sum = sum(sum)) %>%
#   ungroup()

# merge reads with the same TEVs_trans (doesnt do anything here)
data <- data %>%
  group_by(TEVp, TEVs, n_mut_TEVp, n_mut_TEVs, source_TEVp, source_TEVs) %>%
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

# copy to new variable
DMS <- data

# write to file
write.table(data,
  file = paste0(RES_DIR, '/raw_DMS.txt'), # Data_DMS.txt
  quote = F, row.names = F, col.names = T, sep = '\t')

# make line plot of a few variants
for (current_library in 'DMS') {
  print(paste0('Lineplot of library ', current_library))

  # set temp data
  data_lib <- data
  data_lib <- data_lib %>%
    filter(substring(TEVs, 1, 6) == 'ENLYFQ') %>%
    arrange(desc(AUC)) %>%
    mutate(ID_AUC = 1:nrow(.)) %>%
    arrange(desc(sum)) %>%
    mutate(ID_sum = 1:nrow(.)) %>%
    mutate(
      t0 = t1 * 100,
      t1 = t2 * 100,
      t2 = t3 * 100,
      t3 = t4 * 100,
      t4 = t5 * 100) %>%
    select(TEVp, TEVs, ID_AUC, ID_sum, t0, t1, t2, t3, t4, sum, AUC)

  # subset for line plot
  plot_data <- data_lib %>%
    filter(ID_sum <= 100) %>%
    rename(ID = ID_sum)

  plot_data <- plot_data %>%
    select(ID, t0, t1, t2, t3, t4) %>%
    pivot_longer(
      cols = -ID,
      names_to = 'tp',
      values_to = 'fraction') %>%
    mutate(tp = gsub('t', '', tp)) %>%
    mutate(tp = as.numeric(tp))

  # plot empty graph
  plot_melted_empty <- plot_data %>% filter(ID == 960)

  p_line_empty <- ggplot(plot_melted_empty, aes(x = tp, y = fraction, group = ID)) + 
    geom_line(col = 'black') +
    scale_x_continuous('time after induction (h)',
      breaks = seq(0, 5, 1), labels = 0:5, expand = c(0,0)) +
    scale_y_continuous('fraction flipped (%)',
      limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0,0)) +
    theme_SH() +
    theme(legend.position = 'none') +
    coord_cartesian(clip = 'off')

  ggsave(paste0(PLT_DIR, '/01_dotplot_lib', current_library, '_empty.pdf'),
    plot = p_line_empty, width = 3, height = 3, units = c('in'), scale = 1)

  # plot data
  p_line <- ggplot(plot_data, aes(x = tp, y = fraction, group = ID)) + 
    geom_line(color = 'black', alpha = 0.5) +
    scale_x_continuous('time after induction (h)',
      breaks = seq(0, 5, 1), labels = 0:5, expand = c(0,0)) +
    scale_y_continuous('fraction flipped (%)',
      limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0,0)) +
    theme_SH() +
    theme(legend.position = 'none') +
    coord_cartesian(clip = 'off')
  
  ggsave(paste0(PLT_DIR, '/01_dotplot_lib_', current_library, '.pdf'),
    plot = p_line, width = 3, height = 3, units = c('in'), scale = 1)

  # plot specific stuff
  plot_data <- data_lib %>%
    filter(TEVs == 'ENLYFQS') %>%
    filter(nchar(TEVp) <= 5) %>%
    filter(substring(TEVp, 1, 4) == 'C151' | TEVp == 'wt') %>%
    mutate(AA = substring(TEVp, 5, 5))

  plot_data$AA[plot_data$TEVp == 'wt'] <- 'C'

  plot_data <- plot_data %>%
    select(AA, t0, t1, t2, t3, t4) %>%
    pivot_longer(
      cols = -AA,
      names_to = 'tp',
      values_to = 'fraction') %>%
    mutate(tp = gsub('t', '', tp)) %>%
    mutate(tp = as.numeric(tp))

  # plot data
  p <- ggplot(plot_data, aes(x = tp, y = fraction, group = AA, color = AA)) + 
    geom_line() +
    scale_x_continuous('time after induction (h)',
      breaks = seq(0, 5, 1), labels = 0:5, expand = c(0,0)) +
    scale_y_continuous('fraction flipped (%)',
      limits = c(0, 50), breaks = seq(0, 100, 25), expand = c(0,0)) +
    theme_SH() +
    coord_cartesian(clip = 'off')
  
  ggsave(paste0(PLT_DIR, '/01_dotplot_lib_', current_library, '.pdf'),
    plot = p_line, width = 3, height = 3, units = c('in'), scale = 1)

}



##################### MLD VARIANTS #####################
# extract library
data <- raw_MLD %>%
  inner_join(., ref_TEVp_MLD %>% filter(list_TEVp == 'WL'), by = 'BC_TEVp') %>%
  inner_join(., ref_TEVs_controls, by = 'BC_TEVs')

# merge reads with the same TEVs_trans (doesnt do anything here)
data <- data %>%
  group_by(TEVp_ABC, TEVp_DEF, TEVp_XYZ, TEVs, library, source_TEVp, source_TEVs) %>%
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
  file = paste0(RES_DIR, '/raw_MLD.txt'), # Data_MLD.txt
  quote = F, row.names = F, col.names = T, sep = '\t')



##################### DMS HEATMAP #####################
# first filter for single mutants
data <- DMS %>%
  filter(n_mut_TEVp <= 1, n_mut_TEVs <= 1)

# clean
data <- data %>%
  select(TEVp, TEVs, AUC)

# extract WT
data_WT <- data %>% filter(TEVp == 'wt')

codon_AA <- data %>%
  filter(TEVp != 'wt') %>%
  mutate(codon = readr::parse_number(TEVp)) %>%
  mutate(original = substring(TEVp, 1, 1)) %>%
  group_by(original, codon) %>%
  summarize() %>%
  ungroup() %>%
  arrange(codon) %>%
  mutate(TEVp_actual = paste0(original, codon, original)) %>%
  mutate(TEVp = 'wt') %>%
  select(-original, -codon)

data_WT <- data_WT %>%
  inner_join(., codon_AA, by = 'TEVp') %>%
  select(-TEVp) %>%
  rename(TEVp = TEVp_actual)

# add to data
data <- rbind(data, data_WT)
  
# add codon
data <- data %>%
  filter(TEVp != 'wt') %>%
  mutate(codon = readr::parse_number(TEVp)) %>%
  mutate(L = nchar(TEVp)) %>%
  mutate(AA = substring(TEVp, L, L)) %>%
  filter(codon %in% c(1:234))

# plot
p <- data %>% 
  group_by(AA) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  ggplot(., aes(x = AA, y = count)) +
    geom_bar(position = 'dodge', stat = 'identity', fill = 'dodgerblue3', color = 'black') +
    theme_SH() +
    scale_x_discrete('Exchanged TEVp AA') +
    scale_y_continuous('frequency') +
    coord_cartesian(clip = 'off')

ggsave(plot = p, paste0(PLT_DIR, '/DMS_TEVp_AA_distribution.png'),
  width = 6, height = 3)

# filter for WT TEVs
DMS_WT <- data %>%
  filter(TEVs == 'ENLYFQS')

WT_activity <- DMS %>%
  filter(TEVs == 'ENLYFQS') %>%
  filter(TEVp == 'wt') %>%
  pull(AUC)

p <- ggplot(DMS_WT, aes(x = codon, y = AA, fill = AUC)) + 
  geom_tile() +
  scale_fill_distiller('Activity (AUC)', palette = 'YlGnBu') + 
  theme_SH() +
  scale_x_continuous('TEVp position',
    expand = c(0, 0), breaks = seq(0, 240, 10)) + 
  scale_y_discrete('TEVp amino acid',
    expand = c(0, 0)) +
  coord_cartesian(clip = 'off')

ggsave(plot = p, paste0(PLT_DIR, '/DMS_TEVs_WT_heatmap.png'),
  width = 12, height = 3)













###############################################################################
################################ OLD STUFF !!! ################################
###############################################################################






# plot overlap betweem PacBio and DMS in MLDEEP library
PB <- ref_TEVp_MLD %>%
  filter(list == 'WL') %>% filter(!is.na(BC_TEVp)) %>%
  pull(BC_TEVp) %>% unique()

ILL <- raw_MLD %>% filter(!is.na(BC_TEVp)) %>%
  pull(BC_TEVp) %>% unique()

# make overlapping VennDiagram
venn.diagram(
  x = list(PB, ILL),
  category.names = c('PacBio' , 'Illumina'),
  filename = paste0(PLT_DIR, '/BC_TEVp_Venn.png'),
  disable.logging = TRUE,
  OUTPUT = FALSE,
  fill = c('#2D4262', '#D09683'),
  alpha = c(0.5, 0.5),
  lwd = 0
)

# plot MLD overlapping barcodes
COLUMNS <- c(paste0('N', 1:5), paste0('F', 1:5))



# # summarize reads
# data <- data %>%
#   group_by(BC_TEVp, library, TEVs) %>%
#   summarize(
#     N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
#     F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5),
#     sum = sum(sum)) %>%
#   ungroup()

# data$min <- rowMins(as.matrix(data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))]))




























raw_MLD$sum <- rowSums(raw_MLD[, COLUMNS])

data_merged <- 
# rank and plot
data_merged <- data_merged %>%
  arrange(desc(total_reads)) %>%
  mutate(rank = 1:nrow(data_merged))

data_plot <- data_merged %>%
  mutate(cutoff = total_reads > 10)
p <- ggplot(data_plot, aes(x = rank, y = log10(total_reads), color = cutoff)) +
    geom_point(size = 0.2) +
    scale_x_continuous(limits = c(0, 1500), breaks = seq(0, 1500, 100)) +
    theme_SH()

M <- stringdistmatrix(data_plot$BC_TEVp, data_plot$BC_TEVp[1:1100], method =
  'lv', nthread = 10)

min <- rowMins(M)
data_plot$min <- min

p <- ggplot(data_plot %>% filter(min > 2 | min == 0), aes(x = rank, y = log10(total_reads), color = as.factor(min))) +
    geom_point(size = 5, alpha = 0.5) +
    scale_x_continuous(limits = c(0, 1500), breaks = seq(0, 1500, 100)) +
    theme_SH()










































# set working directory
OUT_DIR_HOME='/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/output'
RESULT_DIR='/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/results'
PLOT_DIR='/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/plots'
setwd(RESULT_DIR)

# load libraries
library(tidyverse)
library(stringdist)
library(data.table)
library(matrixStats)
library(ggseqlogo)
library(scico)
library(pals)
library(stringr)
library(spgs)
library(VennDiagram)

# initialize dummy data
data <- data.frame(seq = 'X')

# read data into a list
for (L in c(1:5)) {
  for (R in c(3:4)) {
    for (state in c('N', 'F')) { # FF as F = FALSE
      # set input
      input_dat <- paste0(OUT_DIR_HOME, '/08_L', L, '_R', R, '_', state, '_uniq.txt')
      print(input_dat)
      
      # read data and index
      data_temp <- read.table(input_dat, colClasses = c('character', 'integer'), header = F, sep = '\t')
      names(data_temp) <- c('seq', paste0(state, L, R))

      # merge
      data <- data %>% full_join(., data_temp, by = 'seq')
    }
  }
}

# replace NA with 0
data[is.na(data)] <- 0

# make rowsums
data$sum <- rowSums(data[, -1])

# filter
data_HQ <- data %>% filter(sum >= 5) %>% arrange(desc(sum))

# split into libraries and samples
BC <- read.table('/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/Barcodes.txt',
  header = T, sep = '\t')

data_list <- list()

libs <- unique(BC$lib)

for (idx_lib in libs) {
  # filter BC list for each library
  BC_lib <- BC %>% filter(lib == idx_lib)
  # extract data
  data_temp <- data_HQ[, c('seq',
    paste0('N', BC_lib$L, BC_lib$R),
    paste0('F', BC_lib$L, BC_lib$R)
    )]
  names(data_temp) <- c('seq',
    paste0('N', BC_lib$tp),
    paste0('F', BC_lib$tp)
    )
  data_temp$sum <- rowSums(data_temp[, -1])

  # filter
  data_temp_HQ <- data_temp %>% filter(sum >= 5) %>% arrange(desc(sum))

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
data_all$library[data_all$lib == 2] <- 'MLDEEP'

# find position of SpeI, HindIII and PstI
SpeI <- 'TCTAGT'
HindIII <- 'AAGCTT'
PstI <- 'CTGCAG'

pos_Spe <- as.data.frame(str_locate(data_all$seq, SpeI)) %>%
  rename(start_Spe = start, end_Spe = end)
pos_Hin <- as.data.frame(str_locate(data_all$seq, HindIII)) %>%
  rename(start_Hin = start, end_Hin = end)
pos_Pst <- as.data.frame(str_locate(data_all$seq, PstI)) %>%
  rename(start_Pst = start, end_Pst = end)

data_all <- cbind(data_all, pos_Spe, pos_Hin, pos_Pst)

data_all$BC_TEVp <- substring(data_all$seq, data_all$end_Spe+1, data_all$start_Hin-1)
data_all$BC_TEVs <- substring(data_all$seq, data_all$end_Hin+1, data_all$start_Pst-1)

n1 <- data_all %>% filter(lib == 1) %>% pull(BC_TEVp) %>% nchar() %>%
  table() %>% as.data.frame(stringsAsFactors = F) %>%
  rename(length = '.') %>% mutate(lib = 'DMS', BC = 'TEVp')

n2 <- data_all %>% filter(lib == 2) %>% pull(BC_TEVp) %>% nchar() %>%
  table() %>% as.data.frame(stringsAsFactors = F) %>%
  rename(length = '.') %>% mutate(lib = 'MLDEEP', BC = 'TEVp')

n3 <- data_all %>% filter(lib == 1) %>% pull(BC_TEVs) %>% nchar() %>%
  table() %>% as.data.frame(stringsAsFactors = F) %>%
  rename(length = '.') %>% mutate(lib = 'DMS', BC = 'TEVs')

n4 <- data_all %>% filter(lib == 2) %>% pull(BC_TEVs) %>% nchar() %>%
  table() %>% as.data.frame(stringsAsFactors = F) %>%
  rename(length = '.')%>% mutate(lib = 'MLDEEP', BC = 'TEVs')

n5 <- rbind(n1, n2, n3, n4)
class(n5$length) <- 'integer'

p <- ggplot(n5, aes(x = length, y = Freq/10**6, fill = lib)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  theme_SH() +
  scale_x_continuous('BC length (nt)', limits = c(5, 16),
    breaks = seq(5, 16, 1)) +
  scale_y_continuous('frequency (x 10^6)', limits = c(0, 3)) +
  facet_wrap(~BC, ncol = 1)

ggsave(paste0(PLOT_DIR, '/01_BC_length.png'),
  plot = p, width = 5, height = 3, units = c('in'), scale = 1)

### make data
# clean-up
data <- data_all %>% select(BC_TEVp, BC_TEVs, library,
  paste0('N', 1:5), paste0('F', 1:5))

# merge reads with same BC combinations (ignores reading errors elsewhere)
data <- data %>% group_by(BC_TEVp, BC_TEVs, library) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5)) %>%
  ungroup()

data$sum <- rowSums(data[, c(paste0('N', 1:5), paste0('F', 1:5))])
data$min <- rowMins(as.matrix(data[, c(paste0('N', 1:5))] + data[, c(paste0('F', 1:5))]))

data <- data %>% arrange(desc(sum)) %>%
  filter(!is.na(BC_TEVp)) %>%
  filter(!is.na(BC_TEVs))

# write to file
write.table(data,
  file = '/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/results/00_data.txt',
  quote = F, row.names = F, col.names = T, sep = '\t')

# read table
# data <- read.table('/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/results/00_data.txt',
#   header = T, sep = '\t')

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

p <- data %>% filter(library == 'MLDEEP') %>% group_by(BC_TEVs) %>%
  summarize(n = n()) %>% arrange(desc(n)) %>% mutate(count = 1:nrow(.)) %>%
  mutate(library = 'MLD') %>%
  ggplot(., aes(x = count, y = n, color = library)) +
    geom_point(size = 0.2) +
    theme_SH() +
    scale_x_continuous('sorted TEVs barcodes',
      limits = c(0, 200)) +
    scale_y_continuous('total reads',
      limits = c(0, 100000))

ggsave(paste0(PLOT_DIR, '/01_NGS_TEVs_BC_MLDEEP.png'),
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

raw_DMS <- data %>%
  filter(library == 'DMS') %>%
  full_join(., BC_TEVs_PacBio_filtered, by = 'BC_TEVs') %>%
  full_join(., BC_TEVp_PacBio_filtered, by = 'BC_TEVp')

# merge reads with the same TEVs_trans and TEVp
data_sum <- raw_DMS %>%
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

################ MLDEEP ################
df_MLDEEP <- data %>%
  filter(library == 'MLDEEP') # %>%
  # left_join(., BC_TEVp_PacBio, by = 'BC_TEVp') ### DATA MISSING!!!

# map back TEVs - first step find unique TEVs
MAX_ERROR <- 1

df_correct <- df_MLDEEP %>% filter(BC_TEVs %in% BC_TEVs_manual$BC_TEVs)
df_unmapped <- df_MLDEEP %>% filter(!(BC_TEVs %in% BC_TEVs_manual$BC_TEVs))

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
df_MLDEEP_mapped <- rbind(df_correct, df_mapped, df_unmapped)

p <- df_MLDEEP_mapped %>% group_by(TEVs_trans) %>% summarize(n = n()) %>%
  arrange(desc(n)) %>% mutate(fraction = n / sum(n)) %>%
  mutate(P1prime = substring(TEVs_trans, 7, 7)) %>%
  ggplot(., aes(x = P1prime, y = fraction)) +
  geom_bar(position = 'dodge', stat = 'identity',
    color = 'black', fill = 'dodgerblue3') +
  theme_SH() +
  scale_y_continuous('fraction of all reads',
    limits = c(0, 0.25), expand = c(0, 0))

ggsave(paste0(PLOT_DIR, '/02_MLDEEP_P1_mapped.png'),
  plot = p, width = 4, height = 3, units = c('in'), scale = 1)

p <- df_MLDEEP_mapped %>% group_by(BC_TEVs) %>%
  summarize(n = n()) %>% arrange(desc(n)) %>% mutate(count = 1:nrow(.)) %>%
  mutate(library = 'MLD') %>%
  ggplot(., aes(x = count, y = n, color = library)) +
    geom_point(size = 0.2) +
    theme_SH() +
    scale_x_continuous('sorted TEVs barcodes',
      limits = c(0, 200)) +
    scale_y_continuous('total reads',
      limits = c(0, 125000), breaks = seq(0, 125000, 25000))

ggsave(paste0(PLOT_DIR, '/01_NGS_TEVs_BC_MLDEEP_mapped.png'),
  plot = p, width = 4, height = 3, units = c('in'), scale = 1)


# temp <- df_MLDEEP_mapped %>% group_by(BC_TEVp) %>% summarize(sum = sum(sum)) %>%
#   arrange(desc(sum)) %>% mutate(count = 1:nrow(.))

# p <- ggplot(temp, aes(x = count, y = sum)) +
#   geom_point(size = 0.2, color = 'dodgerblue3') +
#   theme_SH() +
#   scale_x_continuous('sorted TEVp barcodes',
#     limits = c(0, 200000)) +
#   scale_y_continuous('total reads',

# ggsave(paste0(PLOT_DIR, '/01_NGS_TEVp_BC_MLDEEP.png'),
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





















################ MLDEEP INTERNAL CONTROLS ################
# extract controls
df_cont <- data %>%
  filter(library == 'MLDEEP') %>%
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
  file = '/home/hoelsimo/02_Sequencing/NGS/2022_LH_SH_S4_LH/results/01_MLDEEP_internal_controls.txt',
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




raw_DMS <- raw_DMS %>% group_by(mutation_TEVp, trans_TEVs) %>%
  summarize(
    N1 = sum(N1), N2 = sum(N2), N3 = sum(N3), N4 = sum(N4), N5 = sum(N5),
    F1 = sum(F1), F2 = sum(F2), F3 = sum(F3), F4 = sum(F4), F5 = sum(F5))

raw_DMS$min <- rowMins(as.matrix(raw_DMS[, c(paste0('N', 1:5))] + raw_DMS[, c(paste0('F', 1:5))]))
raw_DMS$sum <- rowSums(as.matrix(raw_DMS[, c(paste0('N', 1:5), paste0('F', 1:5))]))

# calculate fraction flipped
for (tp in c(1:5)) {
  raw_DMS[, paste0('t', tp)] <- raw_DMS[, paste0('F', tp)] / 
    (raw_DMS[, paste0('F', tp)] + raw_DMS[, paste0('N', tp)])
}

# filter
data_HQ <- raw_DMS %>% filter(min >= 10) %>% arrange(desc(sum))

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
  write.table(., file = '01_raw_DMS.txt',
    quote = F, row.names = F, col.names = T, sep = '\t')


# check some varinats
data_HQ %>% select(mutation_TEVp, trans_TEVs, min, sum, AUC, t1) %>%
  arrange(desc(AUC)) %>% as.data.frame() %>% head(., 30)


# define AUC function
calculate_auc <- function(dataframe) {
  # copy
  temp <- dataframe

  # calculate fraction flipped
  for (tp in c(1:5)) {
    temp[, paste0('t', tp)] <- temp[, paste0('F', tp)] / 
      (temp[, paste0('F', tp)] + temp[, paste0('N', tp)])
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

  temp$AUC <- round(rowSums(m_auc / 240), 4)

  # add back to original frame
  dataframe <- dataframe %>%
    mutate(AUC = temp$AUC)
  # return
  return(dataframe)
}
