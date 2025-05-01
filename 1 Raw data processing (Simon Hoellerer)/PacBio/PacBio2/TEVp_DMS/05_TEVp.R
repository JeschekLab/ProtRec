# set variables
ROOT_DIR <- '/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio2/TEVp'
RAW_DIR <- paste0(ROOT_DIR, '/raw_data')
TMP_DIR <- paste0(ROOT_DIR, '/tmp')
OUT_DIR <- paste0(ROOT_DIR, '/output')
RES_DIR <- paste0(ROOT_DIR, '/results')
PLT_DIR <- paste0(ROOT_DIR, '/plots')
CCS_DIR <- paste0(RES_DIR, '/process_ccs')

# change directory
setwd(PLT_DIR)

# load libraries
library(tidyverse)
library(ggseqlogo)
library(spgs)
library(bioseq)

# read data
CCS <- read.table(paste0(CCS_DIR, '/processed_ccs.csv'), sep = ',', header = T)
BCs <- read.table(paste0(TMP_DIR, '/03_BCs.txt'),
  sep = '\t', header = T)

# remove @ sign
BCs$ID <- gsub('@', '', BCs$ID)

# rename
CCS <- CCS %>% rename(ID = query_name) %>% select(-target, -TEVp_accuracy)

# merge
df <- CCS %>% inner_join(BCs, by = 'ID') %>% select(-ID) %>%
  rename(mut = TEVp_mutations)

# plot BC length
df$n <- nchar(df$BC)

p <- ggplot(data = df, aes(x = n)) +
  geom_histogram(binwidth = 1, fill = 'dodgerblue3', color = 'black') +
  theme_SH() +
  scale_x_continuous('BC TEVp length') +
  scale_y_continuous('frequency') +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/BC_TEVp_hist.png'), plot = p, width = 3,
  height = 3, units = c('in'), scale = 1)

# plot BCs
temp <- df %>% filter(nchar(BC) == 15)
p <- ggplot() +
  geom_logo(data = temp$BC, method = 'probability', seq_type = 'dna') +
  theme_SH() +
  scale_x_continuous('randomised position in BC TEVp') +
  scale_y_continuous('base fraction', limits = c(0, 1)) +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/BC_TEVp_seq.png'), plot = p, width = 5,
  height = 2, units = c('in'), scale = 1)

# find indels
df$insertion <- grepl('ins', df$mut)
df$deletion <- grepl('del', df$mut)

# translate mutations into actual mutations
maxcols <- df %>% pull(mut) %>% str_count(., pattern = ' ') %>% max() + 1
columns <- paste0('mut', 1:maxcols)

df_separated <- df %>%
  separate(., into = columns, col = mut, remove = F)

# replace with empty string
df_separated[is.na(df_separated)] <- ''

# initiate indel column
df_separated$ins <- ''
df_separated$del <- ''

# cycle through columns, copy indels to end and remove
for (current_column in rev(columns)) {
  # find indels
  with_ins <- grepl('ins', df_separated[, current_column])
  with_del <- grepl('del', df_separated[, current_column])

  # add indel to indel column
  df_separated$ins[with_ins] <-
    paste0(df_separated[, current_column][with_ins], ' ', df_separated$ins[with_ins])
  df_separated$del[with_del] <-
    paste0(df_separated[, current_column][with_del], ' ', df_separated$del[with_del])
  
  # overwride
  df_separated[, current_column][with_ins] <- ''
  df_separated[, current_column][with_del] <- ''
}

# put back mut column
df_separated$mut <- apply(df_separated[, columns], 1, paste, collapse = ' ')

# remove trailing and leading white space and separate again
df_separated$mut <- trimws(df_separated$mut, which = 'both')
df_separated$ins <- trimws(df_separated$ins, which = 'both')
df_separated$del <- trimws(df_separated$del, which = 'both')

# remove double spaces
df_separated <- df_separated %>% 
   mutate(mut = str_replace_all(mut, pattern = '  ', replacement = ' ')) %>% 
   mutate(mut = str_replace_all(mut, pattern = '  ', replacement = ' ')) %>% 
   mutate(mut = str_replace_all(mut, pattern = '  ', replacement = ' ')) %>% 
   mutate(mut = str_replace_all(mut, pattern = '  ', replacement = ' ')) %>% 
   mutate(mut = str_replace_all(mut, pattern = '  ', replacement = ' '))

# remove old columns
df_separated[, columns] <- NULL

# split again
maxcols <- df_separated %>% pull(mut) %>% str_count(., pattern = ' ') %>% max() + 1
columns <- paste0('mut', 1:maxcols)
df_separated <- df_separated %>%
  separate(., into = columns, col = mut, remove = F)

# always remove first letter
for (i in columns) {
  df_separated[, i] <- substring(df_separated[, i], 2, 10)
}

# split in numbers and letters
df_base <- df_separated
for (i in columns) {
  df_base[, i] <- gsub('[[:digit:]]', '', df_base[, i])
}

df_posi <- df_separated
for (i in columns) {
  df_posi[, i] <- as.integer(gsub('[^[:digit:]]', '', df_posi[, i]))
}

# add TEVp0
TEVp <- 'ATGGGCGAAAGCCTGTTTAAAGGTCCGCGTGATTATAATCCGATTAGCAGCAGCATTTGCCATCTGACCAATGAAAGTGATGGTCATACCACCAGTCTGTATGGTATTGGTTTTGGTCCGTTTATTATCACCAATAAACACCTGTTTCGTCGCAATAATGGCACCCTGGTGGTTCAGAGCCTGCATGGTGTTTTTAAAGTTAAAGATACCACCACCCTGCAGCAGCATCTGGTTGATGGTCGTGATATGATTATTATTCGTATGCCGAAAGATTTTCCGCCTTTTCCGCAGAAACTGAAATTTCGTGAACCGCAGCGTGAAGAACGTATTTGTCTGGTTACCACCAATTTTCAGACCAAAAGCATGAGCAGCATGGTTAGCGATACCAGCTGTACCTTTCCGAGCGGTGATGGTATTTTTTGGAAACATTGGATTCAGACCAAAGATGGTCAGTGTGGTAGTCCGCTGGTTAGCACCCGTGATGGTTTTATTGTTGGTATTCATAGCGCCAGCAACTTTACCAATACCAACAACTATTTTACCAGCGTGCCGAAAAACTTCATGGAACTGCTGACCAATCAAGAGGCACAGCAGTGGGTTAGCGGTTGGCGTCTGAATGCAGATAGCGTTCTGTGGGGTGGTCATAAAGTTTTTATGGTGAAACCGGAAGAACCGTTTCAGCCGGTTAAAGAAGCGACCCAGCTGTAA'
TEVp_trans <- seq_translate(dna(TEVp))
df_separated$TEVp <- TEVp

# remove all NA values
df_posi[is.na(df_posi)] <- 0
df_base[is.na(df_base)] <- '0'

# change TEVp sequence
for (i in columns) {
  substring(df_separated$TEVp, df_posi[, i], df_posi[, i]) <-
    df_base[, i]
}

# count
unique(df_separated$TEVp) %>% length() # 40302

# lets translate
df_separated$trans <- seq_translate(dna(df_separated$TEVp))

# count
unique(df_separated$trans) %>% length() # 33232

# find codon number
df_codo <- df_posi
for (i in columns) {
  print(i)
  df_codo[, i] <- ceiling(df_codo[, i] / 3)
}

# find number of unique entries
df_codo$n_mut <- NA
for (i in 1:nrow(df_codo)) {
  df_codo$n_mut[i] <- length(unique(unlist(df_codo[i, columns]))) - 1
}

# count
p <- df_codo %>% group_by(n_mut) %>% summarise(n = n()) %>%
  mutate(frac = n / sum(n) * 100) %>%
  ggplot(data = ., aes(x = n_mut, y = frac)) +
    geom_bar(stat = 'identity', fill = 'dodgerblue3', colour = 'black') +
    scale_x_continuous('number of TEVp AA mutations') +
    scale_y_continuous('fraction (%)',
      limits = c(0, 100), breaks = seq(0, 100, 20)) +
    theme_SH() +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/TEVp_mutations.png'), plot = p, width = 3,
  height = 3, units = c('in'), scale = 1)

# clean up
df_new <- df_codo %>% select(BC, insertion, deletion, ins, del,
  mut, all_of(columns), n_mut)
df_new$trans <- df_separated$trans

# combine into final data frame
df_new$trans <- df_separated$trans

# loop through mutations and extract
for (current_column in columns) {
  print(paste(current_column))
  rows <- df_new[, current_column] != '0'
  old <- substring(TEVp_trans, df_new[rows, current_column], df_new[rows, current_column])
  codon <- df_new[rows, current_column] - 1
  new <- substring(df_new$trans[rows], df_new[rows, current_column], df_new[rows, current_column])
  comb <- paste0(old, codon, new)
  df_new[rows, current_column] <- comb
}

# loop through rows and find unique mutations
df_0 <- df_new %>% filter(n_mut == 0) %>% mutate(TEVp = 'wt')
df_1 <- df_new %>% filter(n_mut == 1) %>% mutate(TEVp = mut1)
df_2 <- df_new %>% filter(n_mut >= 2)

for (i in 1:nrow(df_2)) {
  print(i)
  mutation <- df_2[i, columns] %>% unlist() %>% unique() %>%
    sort() %>% .[-1] %>% paste(., collapse = '_')
  df_2$TEVp[i] <- mutation
}

# make final df
df_final <- rbind(df_0, df_1, df_2) %>%
  select(BC, TEVp, n_mut, insertion, deletion, ins, del, mut) %>%
  rename(BC_TEVp = BC) %>%
  mutate(indel = insertion | deletion) %>%
  mutate(ins = trimws(ins, which = 'both')) %>%
  mutate(del = trimws(del, which = 'both')) %>%
  mutate(which_indel = paste0(ins, ' ', del)) %>%
  select(-ins, -del, -insertion, -deletion) %>%
  mutate(comb = paste0(BC_TEVp, '_', TEVp))

nrow(df_final) # 731714 reads
length(unique(df_final$TEVp)) # 37821 unique TEVp
length(unique(df_final$BC_TEVp)) # 103362 unique BCs

# split into with indel and without indel
df_noind <- df_final %>% filter(!indel)
df_indel <- df_final %>% filter(indel)

# add variants from indel data set if they are also in the non-indel set back
df_indel$alsonoindel <- df_indel$comb %in% df_noind$comb

df_indel$alsonoindel <- df_indel$comb %in% df_noind$comb
df_noind <- rbind(
  df_noind,
  df_indel %>% filter(alsonoindel) %>% select(-alsonoindel)
  )

# remove variants from indel data set if they are also in the non-indel set
df_indel <- df_indel %>% filter(!alsonoindel) %>% select(-alsonoindel)

# check if total appearance is same as total reads
df_noind <- df_noind %>%
  group_by(BC_TEVp, TEVp, n_mut, mut, indel, which_indel) %>%
  mutate(total_app = n()) %>% ungroup() %>%
  group_by(comb) %>% mutate(app= n()) %>% ungroup() %>%
  mutate(match = total_app == app)

sum(!df_noind$match) # 542946 - those variants appear with different indels

df_indel <- df_indel %>%
  group_by(BC_TEVp, TEVp, n_mut, mut, indel, which_indel) %>%
  mutate(total_app = n()) %>% ungroup() %>%
  group_by(comb) %>% mutate(app= n()) %>% ungroup() %>%
  mutate(match = total_app == app)

sum(!df_indel$match) # 4087 - those variants appear with different indels

# add those that deviate into no-indel data frame
df_noind <- rbind(
  df_noind,
  df_indel %>% filter(!match)
  )

df_indel <- df_indel %>% filter(match)

# correct indel column
df_noind <- df_noind %>% mutate(which_indel = '') %>% mutate(indel = FALSE)

# now combine again
df_final <- rbind(df_noind, df_indel)

df_final <- df_final %>%
  group_by(BC_TEVp, TEVp, n_mut, indel, which_indel) %>%
  mutate(total_app = n()) %>%
  ungroup() %>%
  group_by(comb) %>% mutate(app= n()) %>% ungroup() %>%
  mutate(match = total_app == app)

sum(!df_final$match) # 0

# calculate reads again
df_final <- df_final %>% group_by(BC_TEVp, TEVp, n_mut, indel) %>%
  summarize(reads_TEVpxBC = n()) %>%
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


write.table(df_final,
  file = paste0(RES_DIR, '/BC_TEVp_PacBio.txt'),
  quote = F, row.names = F, col.names = T, sep = '\t')

# find unique combinations
w <- df_final %>% group_by(BC_TEVp) %>% summarize(reads_BC_TEVp = n()) %>%
  ungroup() %>% arrange(desc(reads_BC_TEVp))

x <- df_final %>% group_by(BC_TEVp, TEVp, n_mut, indel) %>%
  summarize(reads_TEVpxBC = n()) %>% ungroup() %>%
  arrange(desc(reads_TEVpxBC))

# find number of unique TEVp per BC
y <- x %>% group_by(BC_TEVp) %>%
  summarize(TEVp_per_BC = n()) %>% ungroup() %>%
  arrange(desc(TEVp_per_BC))

# plot number of variants per BC
p <- ggplot(data = y, aes(x = TEVp_per_BC)) +
  geom_histogram(binwidth = 1, fill = 'dodgerblue3', color = 'black') +
  theme_SH() +
  scale_x_continuous('Number of TEVp per BC') +
  scale_y_continuous('frequency') +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/BC_TEVp_hist.png'), plot = p, width = 3,
  height = 3, units = c('in'), scale = 1)

z <- x %>% group_by(TEVp) %>%
  summarize(BCs_per_TEVp = n()) %>% ungroup() %>%
  arrange(desc(BCs_per_TEVp))

df_comb <- x %>%
  inner_join(w, by = 'BC_TEVp') %>%
  inner_join(y, by = 'BC_TEVp') %>%
  inner_join(z, by = 'TEVp')

# calculate fraction
df_comb <- df_comb %>% mutate(fraction_BC_TEVp = reads_TEVpxBC / reads_BC_TEVp)

p <- df_comb %>% group_by(BCs_per_TEVp) %>% summarize(n = n()) %>%
  mutate(n2 = n/BCs_per_TEVp) %>% filter(BCs_per_TEVp <= 19) %>%
  ggplot(data = ., aes(x = BCs_per_TEVp, y = n2)) +
    geom_bar(position = 'dodge', stat = 'identity', fill = 'dodgerblue3', color = 'black') +
    theme_SH() +
    scale_x_continuous('Number of BC per TEVp') +
    scale_y_continuous('frequency') +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/TEVp_BC_hist.png'), plot = p, width = 3,
  height = 3, units = c('in'), scale = 1)

# plot histogram
p <- ggplot(data = df_comb, aes(x = fraction_BC_TEVp)) +
  geom_histogram(binwidth = 0.05, fill = 'dodgerblue3', color = 'black') +
  theme_SH() +
  scale_x_continuous('Ratio of TEVp BC (white/blacklist)') +
  scale_y_continuous('frequency') +
  coord_cartesian(clip = 'off')

ggsave(paste0(PLT_DIR, '/ratio_hist.png'), plot = p, width = 3,
  height = 3, units = c('in'), scale = 1)
