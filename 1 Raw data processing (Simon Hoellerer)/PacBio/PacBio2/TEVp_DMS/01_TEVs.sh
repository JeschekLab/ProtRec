#!/bin/bash

# set variables
ROOT_DIR="/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio2/TEVs"
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
gzip -dc "${RAW_DIR}/m64156_221027_051256.hifi_reads.fastq.gz" \
  > "${RAW_DIR}/CCS.fastq"

# split into 4 files for each lane
awk 'NR%4==1' "${RAW_DIR}/CCS.fastq" > "${TMP_DIR}/CCS_01.txt"
awk 'NR%4==2' "${RAW_DIR}/CCS.fastq" > "${TMP_DIR}/CCS_02.txt"
awk 'NR%4==3' "${RAW_DIR}/CCS.fastq" > "${TMP_DIR}/CCS_03.txt"
awk 'NR%4==0' "${RAW_DIR}/CCS.fastq" > "${TMP_DIR}/CCS_04.txt"

# extract TEVs data based on constant TEVp barcode
agrep -3 "CGTAAAACTAGTAGCAGCCTACCGAAGCTT" "${TMP_DIR}/CCS_02.txt" \
  > "${TMP_DIR}/data_TEVs.txt"

# add reverse complement
agrep -3 "AAGCTTCGGTAGGCTGCTACTAGTTTTACG" "${TMP_DIR}/CCS_02.txt" \
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
