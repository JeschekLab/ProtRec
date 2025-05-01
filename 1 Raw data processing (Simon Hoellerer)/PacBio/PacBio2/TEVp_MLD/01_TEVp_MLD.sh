#!/bin/bash

# set variables
ROOT_DIR="/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio2/TEVp_MLD"
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
gzip -dc "${RAW_DIR}/m64156_221027_051256.hifi_reads.fastq.gz" \
  > "${RAW_DIR}/CCS.fastq"

# split into 4 files for each lane
awk 'NR%4==1' "${RAW_DIR}/CCS.fastq" > "${TMP_DIR}/01_CCS_01.txt" 
awk 'NR%4==2' "${RAW_DIR}/CCS.fastq" > "${TMP_DIR}/01_CCS_02.txt"
awk 'NR%4==3' "${RAW_DIR}/CCS.fastq" > "${TMP_DIR}/01_CCS_03.txt"
awk 'NR%4==0' "${RAW_DIR}/CCS.fastq" > "${TMP_DIR}/01_CCS_04.txt"

# add reverse
cat "${TMP_DIR}/01_CCS_02.txt" | rev | tr ACGT TGCA > "${TMP_DIR}/01_CCS_02_rev.txt" 
cat "${TMP_DIR}/01_CCS_04.txt" | rev > "${TMP_DIR}/01_CCS_04_rev.txt"

# combine forward and reverse
cat "${TMP_DIR}/01_CCS_02.txt" "${TMP_DIR}/01_CCS_02_rev.txt" \
  > "${TMP_DIR}/02_CCS.txt"

# extract MLD reads with constant BC_TEVs (and only 2 MM)
agrep -2 -V0 "CTGCAGACATCGAAGCTT" \
  "${TMP_DIR}/02_CCS.txt" > "${TMP_DIR}/03_CCS_MLD.txt"

### find constant sequences
# find upstream of barcode
cat "${TMP_DIR}/03_CCS_MLD.txt" \
  | tre --max-errors=1 -n --show-position "CACTGCAGACATCGAAGCTT" \
  | tr ':' '\t' | tr '-' '\t' \
  | awk 'BEGIN {OFS = "\t"} {print $1, $2, $3}' \
  > "${OUT_DIR}/BC_L.txt" & 

# find downstream of barcode
cat "${TMP_DIR}/03_CCS_MLD.txt" \
  | tre --max-errors=1 -n --show-position "ACTAGTTTTACGGCTAGCTC" \
  | tr ':' '\t' | tr '-' '\t' \
  | awk 'BEGIN {OFS = "\t"} {print $1, $2, $3}' \
  > "${OUT_DIR}/BC_R.txt" &

# find upstream of ABC
cat "${TMP_DIR}/03_CCS_MLD.txt" \
  | tre --max-errors=1 -n --show-position "ATGAAAGTGATGGTCATACC" \
  | tr ':' '\t' | tr '-' '\t' \
  | awk 'BEGIN {OFS = "\t"} {print $1, $2, $3}' \
  > "${OUT_DIR}/ABC_L.txt" &

# find downstream of ABC
cat "${TMP_DIR}/03_CCS_MLD.txt" \
  | tre --max-errors=1 -n --show-position "TATGGTATTGGTTTTGGTCC" \
  | tr ':' '\t' | tr '-' '\t' \
  | awk 'BEGIN {OFS = "\t"} {print $1, $2, $3}' \
  > "${OUT_DIR}/ABC_R.txt" &

# find upstream of DEF
cat "${TMP_DIR}/03_CCS_MLD.txt" \
  | tre --max-errors=1 -n --show-position "AACATTGGATTCAGACCAAA" \
  | tr ':' '\t' | tr '-' '\t' \
  | awk 'BEGIN {OFS = "\t"} {print $1, $2, $3}' \
  > "${OUT_DIR}/DEF_L.txt" &

# find downstream of DEF
cat "${TMP_DIR}/03_CCS_MLD.txt" \
  | tre --max-errors=1 -n --show-position "TGTGGTAGTCCGCTGGTTAG" \
  | tr ':' '\t' | tr '-' '\t' \
  | awk 'BEGIN {OFS = "\t"} {print $1, $2, $3}' \
  > "${OUT_DIR}/DEF_R.txt" &

# find upstream of XYZ
cat "${TMP_DIR}/03_CCS_MLD.txt" \
  | tre --max-errors=1 -n --show-position "TGTGGGGTGGTCATAAAGTT" \
  | tr ':' '\t' | tr '-' '\t' \
  | awk 'BEGIN {OFS = "\t"} {print $1, $2, $3}' \
  > "${OUT_DIR}/XYZ_L.txt" &

# find downstream of XYZ
cat "${TMP_DIR}/03_CCS_MLD.txt" \
  | tre --max-errors=1 -n --show-position "AAACCGGAAGAACCGTTTCA" \
  | tr ':' '\t' | tr '-' '\t' \
  | awk 'BEGIN {OFS = "\t"} {print $1, $2, $3}' \
  > "${OUT_DIR}/XYZ_R.txt" &
wait
