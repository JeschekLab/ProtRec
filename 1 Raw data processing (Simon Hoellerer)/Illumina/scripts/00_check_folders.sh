#!/bin/bash

# import config
source "${CONFIG}"

# create folders and check integrity
for current_dir in ${RAW_DIR} ${OUT_DIR} ${TMP_DIR} ${LOG_DIR} ${RES_DIR} ${PLT_DIR}; do
  # create folders and check integrity
  if [ ! -d "${current_dir}" ]; then
    mkdir -p "${current_dir}"
    if [ ! -d "${current_dir}" ]; then
      printf "ERROR: ${current_dir} does not exist and cannot be created!\n"
      printf "Directory = ${current_dir}\n"
      printf "Exiting the script.\n"
      exit
    fi
    printf "${current_dir} successfully created. Continuing ...\n"
  else
    printf "${current_dir} already exists. Continuing ...\n"
  fi
done
