#!/bin/bash

# set variables
ROOT_DIR="/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/PacBio2/TEVp"
RAW_DIR="${ROOT_DIR}/raw_data"
TMP_DIR="${ROOT_DIR}/tmp"
DAT_DIR="${ROOT_DIR}/data"
OUT_DIR="${ROOT_DIR}/output"
RES_DIR="${ROOT_DIR}/results"
PLT_DIR="${ROOT_DIR}/plots"

# pack gz files
gzip -c "${TMP_DIR}/03_CCS_cropped.fastq" > "${DAT_DIR}/CCS_cropped.fastq.gz"
