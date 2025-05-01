# set variables
ROOT_DIR <- "/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio2/TEVp"
RAW_DIR <- paste0(ROOT_DIR, "/raw_data")
TMP_DIR <- paste0(ROOT_DIR, "/tmp")
OUT_DIR <- paste0(ROOT_DIR, "/output")
RES_DIR <- paste0(ROOT_DIR, "/results")
PLT_DIR <- paste0(ROOT_DIR, "/plots")

# change directory
setwd(TMP_DIR)

# load libraries
library(tidyverse)

# read data
a1 <- read.table("01_CCS_01.txt", sep = "\t", header = F,
  colClasses = "character")
a3 <- read.table("01_CCS_03.txt", sep = "\t", header = F,
  colClasses = "character")
a4_fwd <- read.table("01_CCS_04.txt", sep = "\t", header = F,
  colClasses = "character", quote = "", comment.char = "")
a4_rev <- read.table("01_CCS_04_rev.txt", sep = "\t", header = F,
  colClasses = "character", quote = "", comment.char = "")

df_fwd <- cbind(a1, a3, a4_fwd)
df_rev <- cbind(a1, a3, a4_rev)

names(df_fwd) <- c("V1", "V3", "V4")
names(df_rev) <- c("V1", "V3", "V4")

df_fwd$line <- 1:nrow(df_fwd)
df_rev$line <- 1:nrow(df_rev)

b1 <- read.table("02_term_fwd.txt", sep = "\t", header = F,
  colClasses = c(rep("integer", 3), "character"))
b2 <- read.table("02_spac_fwd.txt", sep = "\t", header = F,
  colClasses = c(rep("integer", 3), "character"))
b3 <- read.table("02_term_rev.txt", sep = "\t", header = F,
  colClasses = c(rep("integer", 3), "character"))
b4 <- read.table("02_spac_rev.txt", sep = "\t", header = F,
  colClasses = c(rep("integer", 3), "character"))

b5 <- read.table("02_BCl_fwd.txt", sep = "\t", header = F,
  colClasses = c(rep("integer", 3), "character"))
b6 <- read.table("02_BCr_fwd.txt", sep = "\t", header = F,
  colClasses = c(rep("integer", 3), "character"))
b7 <- read.table("02_BCl_rev.txt", sep = "\t", header = F,
  colClasses = c(rep("integer", 3), "character"))
b8 <- read.table("02_BCr_rev.txt", sep = "\t", header = F,
  colClasses = c(rep("integer", 3), "character"))

names(b1) <- c("line", "term_fwd_start", "term_fwd_end", "seq")
names(b2) <- c("line", "spac_fwd_start", "spac_fwd_end", "seq")
names(b3) <- c("line", "term_rev_start", "term_rev_end", "seq")
names(b4) <- c("line", "spac_rev_start", "spac_rev_end", "seq")
names(b5) <- c("line", "BCl_fwd_start", "BCl_fwd_end", "seq")
names(b6) <- c("line", "BCr_fwd_start", "BCr_fwd_end", "seq")
names(b7) <- c("line", "BCl_rev_start", "BCl_rev_end", "seq")
names(b8) <- c("line", "BCr_rev_start", "BCr_rev_end", "seq")

# use all data, that has all 4 fwd BCl, BCr, spac and terminator
c1 <- df_fwd %>%
  inner_join(b1, by = "line") %>%
  inner_join(b2, by = "line") %>%
  inner_join(b5, by = "line") %>%
  inner_join(b6, by = "line")

c2 <- df_rev %>%
  inner_join(b3, by = "line") %>%
  inner_join(b4, by = "line") %>%
  inner_join(b7, by = "line") %>%
  inner_join(b8, by = "line")

c1 <- c1 %>%
  select(V1, seq.x, V3, V4, BCl_fwd_end, BCr_fwd_start, term_fwd_start, spac_fwd_end) %>%
  rename(V2 = seq.x, BC_start = BCl_fwd_end, BC_end = BCr_fwd_start,
    term_start = term_fwd_start, spac_end = spac_fwd_end)
c2 <- c2 %>%
  select(V1, seq.x, V3, V4, BCl_rev_end, BCr_rev_start, term_rev_start, spac_rev_end) %>%
  rename(V2 = seq.x,
    BC_start = BCl_rev_end, BC_end = BCr_rev_start,
    term_start = term_rev_start, spac_end = spac_rev_end)

# combine both directions
d <- rbind(c1, c2) # %>% select(-BC_end, -BC_start)

# extract BC sequences
d$BC <- substring(d$V2, d$BC_start+1, d$BC_end)

# extract sequences and quality
d$V2_new <- substring(d$V2, d$term_start+1, d$spac_end)
d$V4_new <- substring(d$V4, d$term_start+1, d$spac_end)

# remove reads with strange BCs
e <- d %>%
  filter(nchar(BC) < 20) %>% filter(nchar(BC) > 10) %>%
  filter(nchar(V2_new) < 2000) %>% filter(nchar(V2_new) > 500) %>%
  filter(nchar(V4_new) < 2000) %>% filter(nchar(V4_new) > 500) %>%
  select(-BC_start, -BC_end, -term_start, -spac_end)

# make fastq again
n <- nrow(e) * 4
f <- data.frame(line = 1:n)
f$line[seq(1, n, 4)] <- e$V1
f$line[seq(2, n, 4)] <- e$V2_new
f$line[seq(3, n, 4)] <- e$V3
f$line[seq(4, n, 4)] <- e$V4_new

# write to file
write.table(f, file = "03_CCS_cropped.fastq", sep = "\t", quote = F,
  col.names = F, row.names = F)

# also extract BCs
g <- e %>% select(V1, BC) %>% rename(ID = V1)

# write to file
write.table(g, file = "03_BCs.txt", sep = "\t", quote = F,
  col.names = T, row.names = F)
