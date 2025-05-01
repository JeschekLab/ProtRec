#!/bin/bash

# set variables
ROOT_DIR="/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio1/TEVs_DMS"
RAW_DIR="${ROOT_DIR}/raw_data"
TMP_DIR="${ROOT_DIR}/tmp"
OUT_DIR="${ROOT_DIR}/output"
RES_DIR="${ROOT_DIR}/results"
PLT_DIR="${ROOT_DIR}/plots"

# create folder
mkdir ${TMP_DIR}
mkdir ${OUT_DIR}
mkdir ${RES_DIR}
mkdir ${PLT_DIR}

# change directory
cd ${RAW_DIR}

# extract file
gzip -dc "${RAW_DIR}/Project_1758_analysis-unknown-814-m54273_211130_005735.Q20.fastq.gz" \
  > "${RAW_DIR}/CCS.fastq"

# split into all 4 lines
awk 'NR%4==0' "CCS.fastq" \
  > "${TMP_DIR}/data_0.txt"
awk 'NR%4==1' "CCS.fastq" \
  > "${TMP_DIR}/data_1.txt"
awk 'NR%4==2' "CCS.fastq" \
  > "${TMP_DIR}/data_2.txt"
awk 'NR%4==3' "CCS.fastq" \
  > "${TMP_DIR}/data_3.txt"

### CHANGE ###
# # split into 4 files for each lane
# awk 'NR%4==1' "${RAW_DIR}/CCS.fastq" > "CCS_01.txt"
# awk 'NR%4==2' "${RAW_DIR}/CCS.fastq" > "CCS_02.txt"
# awk 'NR%4==3' "${RAW_DIR}/CCS.fastq" > "CCS_03.txt"
# awk 'NR%4==0' "${RAW_DIR}/CCS.fastq" > "CCS_04.txt"

# extract TEVs data based on 50 bp from MBP-tag
agrep -3 "CGTAAAACTAGTAGCAGCCTACCGAAGCTT" "${TMP_DIR}/data_2.txt" \
  > "${TMP_DIR}/data_TEVs.txt"

# add reverse complement
agrep -3 "AAGCTTCGGTAGGCTGCTACTAGTTTTACG" "${TMP_DIR}/data_2.txt" \
  | rev | tr ACGT TGCA \
  >> "${TMP_DIR}/data_TEVs.txt"

# find position of TEVs and BC
cat "${TMP_DIR}/data_TEVs.txt" \
  | tre --max-errors=3 -n --show-position "GCTGGCTCCGCTGCTGGTTCTGGCTCGGGC" \
  | tr ':' '\t' | tr '-' '\t' > "${OUT_DIR}/TEVs_left.txt" & 
cat "${TMP_DIR}/data_TEVs.txt" \
  | tre --max-errors=3 -n --show-position "GCTAGCGCATGCGGTGGCAGCGGAGGTTCG" \
  | tr ':' '\t' | tr '-' '\t' > "${OUT_DIR}/TEVs_right.txt" & 
cat "${TMP_DIR}/data_TEVs.txt" \
  | tre --max-errors=3 -n --show-position "CGTAAAACTAGTAGCAGCCTACCGAAGCTT" \
  | tr ':' '\t' | tr '-' '\t' > "${OUT_DIR}/BC_left.txt" & 
cat "${TMP_DIR}/data_TEVs.txt" \
  | tre --max-errors=3 -n --show-position "CTGCAGTGAAAGGATTAGGGGAACGCAGAG" \
  | tr ':' '\t' | tr '-' '\t' > "${OUT_DIR}/BC_right.txt" & 
wait
