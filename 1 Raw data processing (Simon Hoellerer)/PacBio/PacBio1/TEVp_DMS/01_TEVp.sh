#!/bin/bash

# set variables
ROOT_DIR="/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio1/TEVp_DMS"
RAW_DIR="${ROOT_DIR}/raw_data"
TMP_DIR="${ROOT_DIR}/tmp"
DAT_DIR="${ROOT_DIR}/data"
OUT_DIR="${ROOT_DIR}/output"
RES_DIR="${ROOT_DIR}/results"
PLT_DIR="${ROOT_DIR}/plots"

# create folder
mkdir ${TMP_DIR}
mkdir ${OUT_DIR}
mkdir ${RES_DIR}
mkdir ${PLT_DIR}
mkdir ${DAT_DIR}

# change directory
cd ${RAW_DIR}

# extract file
gzip -dc "${RAW_DIR}/Project_1758_analysis-unknown-814-m54273_211130_005735.Q20.fastq.gz" \
  > "${RAW_DIR}/CCS.fastq"

# split into 4 files for each lane
awk 'NR%4==1' "${RAW_DIR}/CCS.fastq" > "${TMP_DIR}/01_CCS_01.txt" 
awk 'NR%4==2' "${RAW_DIR}/CCS.fastq" > "${TMP_DIR}/01_CCS_02.txt"
awk 'NR%4==3' "${RAW_DIR}/CCS.fastq" > "${TMP_DIR}/01_CCS_03.txt"
awk 'NR%4==0' "${RAW_DIR}/CCS.fastq" > "${TMP_DIR}/01_CCS_04.txt"

# add reverse
cat "${TMP_DIR}/01_CCS_02.txt" | rev | tr ACGT TGCA > "${TMP_DIR}/01_CCS_02_rev.txt" 
cat "${TMP_DIR}/01_CCS_04.txt" | rev > "${TMP_DIR}/01_CCS_04_rev.txt" 

# find forward terminator
cat "${TMP_DIR}/01_CCS_02.txt" \
  | tre --max-errors=5 -n --show-position "TCTAGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAA" \
  | tr ':' '\t' | tr '-' '\t' > "${TMP_DIR}/02_term_fwd.txt" & 

# find forward RBS
cat "${TMP_DIR}/01_CCS_02.txt" \
  | tre --max-errors=5 -n --show-position "TAGACCGGTGCTAGCATAATACCTAGGGCTGAGCTAGCCGTAAAACTAGT" \
  | tr ':' '\t' | tr '-' '\t' > "${TMP_DIR}/02_spac_fwd.txt" &

# find reverse terminator
cat "${TMP_DIR}/01_CCS_02_rev.txt" \
  | tre --max-errors=5 -n --show-position "TCTAGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAA" \
  | tr ':' '\t' | tr '-' '\t' > "${TMP_DIR}/02_term_rev.txt" &

# find reverse RBS
cat "${TMP_DIR}/01_CCS_02_rev.txt" \
  | tre --max-errors=5 -n --show-position "ACAAACAATCCCTCAAAGATCGGGCCCTAAAGCTCGCAAGATTTAACCAC" \
  | tr ':' '\t' | tr '-' '\t' > "${TMP_DIR}/02_spac_rev.txt" &

# find forward upstream of BC
cat "${TMP_DIR}/01_CCS_02.txt" \
  | tre --max-errors=3 -n --show-position "ACCTAGGGCTGAGCTAGCCGTAAAACTAGT" \
  | tr ':' '\t' | tr '-' '\t' > "${TMP_DIR}/02_BCl_fwd.txt" &

# find forward downstream of BC
cat "${TMP_DIR}/01_CCS_02.txt" \
  | tre --max-errors=3 -n --show-position "AAGCTTACGTCGCTGACTGCAGTGAAAGGA" \
  | tr ':' '\t' | tr '-' '\t' > "${TMP_DIR}/02_BCr_fwd.txt" &

# find reverse upstream of BC
cat "${TMP_DIR}/01_CCS_02_rev.txt" \
  | tre --max-errors=3 -n --show-position "ACCTAGGGCTGAGCTAGCCGTAAAACTAGT" \
  | tr ':' '\t' | tr '-' '\t' > "${TMP_DIR}/02_BCl_rev.txt" &

# find reverse downstream of BC
cat "${TMP_DIR}/01_CCS_02_rev.txt" \
  | tre --max-errors=3 -n --show-position "AAGCTTACGTCGCTGACTGCAGTGAAAGGA" \
  | tr ':' '\t' | tr '-' '\t' > "${TMP_DIR}/02_BCr_rev.txt" &

wait
