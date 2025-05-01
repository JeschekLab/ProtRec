# import variables
ROOT_DIR <- '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/uASPIre_screen4_5/Illumina_PacBio_combined'
PLT_DIR <- paste0(ROOT_DIR, '/plots')
RES_DIR <- paste0(ROOT_DIR, '/results')
Illumina_DIR <- '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/uASPIre_screen4_5/2022_LH_SH_S4_combined'
PacBio_DIR <- '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio_combined'

# change directory
setwd(ROOT_DIR)

# load libraries
library(Hmisc)
library(tidyverse)
library(stringi)
library(cowplot)
library(scico)

# check positions
# ABC: 30, 31, 32
# DEF: 148, 149, 150
# XYZ: 217, 218, 219 

# read data
data <- read.table(
  file = paste0(RES_DIR, '/Data_DMS.txt'),
  header = T, sep = '\t')

# only keep main columns
data <- data %>%
  filter(n_mut_TEVp == 1) %>%
  filter(n_mut_TEVs == 1) %>%
  select(TEVp, TEVs, sum, min, AUC)

# find out which TEVs is mutated
data <- data %>%
  mutate(TEVs_pos = NA)

# loop through all 7 TEVs positions and find mismatch with wildtype
TEVs <- 'ENLYFQS'
for (current_position in 1:7) {
  ref <- substring(TEVs, current_position, current_position)
  df <- substring(data$TEVs, current_position, current_position)
  match <- !(ref == df)
  data$TEVs_pos[match] <- current_position
}

data <- data %>%
  mutate(TEVs_AA = substring(TEVs, TEVs_pos, TEVs_pos))

# extract TEVp mutation
data <- data %>%
  mutate(TEVp_short = substring(TEVp, 2, 5)) %>%
  mutate(TEVp_rev = stri_reverse(TEVp_short)) %>%
  mutate(TEVp_pos = as.integer(stri_reverse(substring(TEVp_rev, 2, 5)))) %>%
  mutate(TEVp_AA = substring(TEVp_rev, 1, 1)) %>%
  select(-TEVp_rev, -TEVp_short)

# calculate average activity for each TEVp
data_avg_TEVp <- data %>%
  group_by(TEVp) %>%
  summarize(
    mean_TEVp = mean(AUC),
    SD_TEVp = sd(AUC)) %>%
  ungroup()

# calculate average activity for each position in TEVp
data_avg_TEVp_pos <- data %>%
  group_by(TEVp_pos) %>%
  summarize(
    mean_TEVp_pos = mean(AUC),
    SD_TEVp_pos = sd(AUC)) %>%
  ungroup()

# calculate average activity for each AA in TEVp
data_avg_TEVp_AA <- data %>%
  group_by(TEVp_AA) %>%
  summarize(
    mean_TEVp_AA = mean(AUC),
    SD_TEVp_AA = sd(AUC)) %>%
  ungroup()

# calculate average activity for each TEVs
data_avg_TEVs <- data %>%
  group_by(TEVs) %>%
  summarize(
    mean_TEVs = mean(AUC),
    SD_TEVs = sd(AUC)) %>%
  ungroup()

# calculate average activity for each position in TEVs
data_avg_TEVs_pos <- data %>%
  group_by(TEVs_pos) %>%
  summarize(
    mean_TEVs_pos = mean(AUC),
    SD_TEVs_pos = sd(AUC)) %>%
  ungroup()

# calculate average activity for each AA in TEVs
data_avg_TEVs_AA <- data %>%
  group_by(TEVs_AA) %>%
  summarize(
    mean_TEVs_AA = mean(AUC),
    SD_TEVs_AA = sd(AUC)) %>%
  ungroup()

# merge
data <- data %>%
  inner_join(., data_avg_TEVp, by = 'TEVp') %>%
  inner_join(., data_avg_TEVp_pos, by = 'TEVp_pos') %>%
  inner_join(., data_avg_TEVp_AA, by = 'TEVp_AA') %>%
  inner_join(., data_avg_TEVs, by = 'TEVs') %>%
  inner_join(., data_avg_TEVs_pos, by = 'TEVs_pos') %>%
  inner_join(., data_avg_TEVs_AA, by = 'TEVs_AA')

### check if ABC, DEF, XYZ would be selected
data_lib <- data %>%
  filter(TEVs_AA != '*', TEVp_AA != '*') %>%
  filter(TEVp_pos >= 1) %>%
  filter(TEVs_pos == 7) %>%
  group_by(TEVp_pos, TEVp_AA) %>%
  summarize(mean_P1p = mean(AUC, na.rm = T)) %>%
  ungroup() %>%
  # group_by(TEVp_AA) %>%
  mutate(rank = rank(-mean_P1p)) %>%
  ungroup() %>%
  # mutate(rank2 = ifelse(rank > 20, 30, rank)) %>%
  mutate(top20 = rank <= 100)

p <- data %>%
  filter(TEVs_AA != '*', TEVp_AA != '*') %>%
  filter(TEVp_pos >= 1) %>%
  filter(TEVs_pos == 7) %>%
  filter(TEVs_AA == 'P') %>%
  ggplot(., aes(x = TEVp_pos, y = TEVp_AA, fill = AUC)) +
    geom_tile() +
    theme_SH() +
    scale_fill_scico(palette = 'batlow') +
    # scale_fill_manual(values = c('black', 'dodgerblue3')) +
    scale_x_continuous(
      name = 'TEVp position',
      expand = c(0, 0),
      limits = c(0, 240),
      breaks = seq(0, 240, 10)) +
    scale_y_discrete(name = 'TEVp amino acid') +
    coord_cartesian(clip = 'off')

save_plot(
  filename = './plots/DMS_P1p_P.pdf',
  base_width = 20,
  p)

# which TEVp position changes the activity in which TEVs position
data_avg_TEVp_TEVs <- data %>%
  filter(TEVs_AA != '*', TEVp_AA != '*') %>%
  filter(TEVp_pos >= 1) %>%
  group_by(TEVp_pos, TEVs_pos) %>%
  summarize(
    mean_TEVp_TEVs = mean(AUC),
    SD_TEVp_TEVs = sd(AUC)) %>%
  ungroup() %>%
  inner_join(., data_avg_TEVp_pos, by = 'TEVp_pos') %>%
  inner_join(., data_avg_TEVs_pos, by = 'TEVs_pos') %>%
  select(-SD_TEVp_pos, -SD_TEVs_pos) %>%
  # normalize to each TEVs average position
  mutate(norm_TEVs = mean_TEVp_TEVs / mean_TEVs_pos) %>%
  mutate(norm_TEVp = mean_TEVp_TEVs / mean_TEVp_pos) %>%
  mutate(norm_TEVp_TEVs = mean_TEVp_TEVs / mean_TEVp_pos / mean_TEVs_pos) %>%
  group_by(TEVs_pos) %>%
  mutate(min_TEVs = min(norm_TEVp_TEVs)) %>%
  ungroup() %>%
  mutate(rel = norm_TEVp_TEVs - min_TEVs) %>%
  group_by(TEVs_pos) %>%
  mutate(max_TEVs = max(rel)) %>%
  ungroup() %>%
  mutate(rel2 = rel / max_TEVs)

# make ranks
data_temp <- data_avg_TEVp_TEVs %>%
  group_by(TEVs_pos) %>%
  mutate(rank_SD = rank(-SD_TEVp_TEVs)) %>%
  ungroup() %>%
  mutate(rank_SD2 = ifelse(rank_SD > 20, 30, rank_SD))

p <- data_temp %>%
  ggplot(., aes(x = TEVp_pos, y = TEVs_pos, fill = rank_SD2)) +
    geom_tile() +
    theme_SH() +
    scale_fill_scico(palette = 'batlow') +
    # scale_fill_distiller('Activity (AUC)', palette = 'YlGnBu') + 
    scale_x_continuous(
      name = 'TEVp position',
      expand = c(0, 0),
      limits = c(0, 240),
      breaks = seq(0, 240, 10)) +
    scale_y_continuous(
      name = 'TEVs position',
      breaks = 1:7,
      labels = c(paste0('p', 6:1), "p1'")) +
    coord_cartesian(clip = 'off')



# scale_color_npg()


# plot average TEVs pos
p <- data %>%
  filter(TEVs_AA != '*', TEVp_AA != '*') %>%
  ggplot(., aes(x = TEVs_pos, y = AUC, group = TEVs_pos)) +
    geom_violin(fill = 'dodgerblue3', scale = 'width') +
    stat_summary(
      fun = mean,
      geom = 'point',
      color = 'black') +
    theme_SH() +
    scale_x_continuous(
      name = 'TEVs position',
      breaks = 1:7,
      labels = c(paste0('p', 6:1), "p1'")) +
    scale_y_continuous(
      name = 'IFP',
      expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 0.20))

save_plot(
  filename = 'TEVs_pos.pdf',
  p)

# plot average TEVs AA
p <- data %>%
  filter(TEVs_AA != '*', TEVp_AA != '*') %>%
  ggplot(., aes(x = TEVs_AA, y = AUC, group = TEVs_AA)) +
    geom_violin(fill = 'dodgerblue3', scale = 'width') +
    stat_summary(
      fun = mean,
      geom = 'point',
      color = 'black') +
    theme_SH() +
    scale_x_discrete(name = 'TEVs amino acid') +
    scale_y_continuous(
      name = 'IFP',
      expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 0.20))

save_plot(
  filename = 'TEVs_AA.pdf',
  p)

# plot average TEVp pos
p <- data %>%
  filter(TEVs_AA != '*', TEVp_AA != '*') %>%
  ggplot(., aes(x = TEVp_pos, y = AUC, group = TEVp_pos)) +
    geom_violin(fill = 'dodgerblue3', scale = 'width') +
    stat_summary(
      fun = mean,
      geom = 'point',
      color = 'black') +
    theme_SH() +
    scale_x_continuous(
      name = 'TEVp position') +
    scale_y_continuous(
      name = 'IFP',
      expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 0.20))

save_plot(
  filename = 'TEVp_pos.pdf',
  p)

# plot average TEVp AA
p <- data %>%
  filter(TEVs_AA != '*', TEVp_AA != '*') %>%
  ggplot(., aes(x = TEVp_AA, y = AUC, group = TEVp_AA)) +
    geom_violin(fill = 'dodgerblue3', scale = 'width') +
    stat_summary(
      fun = mean,
      geom = 'point',
      color = 'black') +
    theme_SH() +
    scale_x_discrete(name = 'TEVp amino acid') +
    scale_y_continuous(
      name = 'IFP',
      expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 0.20))

save_plot(
  filename = 'TEVp_AA.pdf',
  p)












# plot average TEVs AA
p <- data_avg_TEVs_AA %>%
  filter(TEVs_AA != '*') %>%
  ggplot(., aes(x = TEVs_AA, y = mean_TEVs_AA)) +
    geom_bar(
      position = 'dodge',
      stat = 'identity',
      color = 'black',
      fill = 'dodgerblue3') +
    geom_errorbar(
      aes(ymin = mean_TEVs_AA-SD_TEVs_AA, ymax = mean_TEVs_AA+SD_TEVs_AA),
      width = 0.2) +
    theme_SH() +
    scale_x_discrete(name = 'TEVs amino acid') +
    scale_y_continuous(
      name = 'average AA activity',
      limits = c(0, 0.125),
      expand = c(0, 0)) +
    coord_cartesian(clip = 'off')

save_plot(
  filename = 'TEVs_AA.pdf',
  p)

# plot average TEVs position 
p <- data_avg_TEVs_pos %>%
  ggplot(., aes(x = TEVs_pos, y = mean_TEVs_pos)) +
    geom_bar(
      position = 'dodge',
      stat = 'identity',
      color = 'black',
      fill = 'dodgerblue3') +
    geom_errorbar(
      aes(ymin = mean_TEVs_pos-SD_TEVs_pos, ymax = mean_TEVs_pos+SD_TEVs_pos),
      width = 0.2) +
    theme_SH() +
    scale_x_continuous(
      name = 'TEVs position',
      breaks = 1:7) +
    scale_y_continuous(
      name = 'average AA activity',
      limits = c(0, 0.15),
      expand = c(0, 0)) +
    coord_cartesian(clip = 'off')

save_plot(
  filename = 'TEVs_AA.pdf',
  p)
