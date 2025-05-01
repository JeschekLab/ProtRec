#!/bin/bash

###############################################################################
# uASPIRE                                                                     #
# Bash executable for NGS raw data analysis                                   #
# Simon Höllerer, ETH Zürich, D-BSSE, BPL, Basel, Switzerland                 #
# E-mail: simon.hoellerer@bsse.ethz.ch                                        #
# Authors: Simon Höllerer                                                     #
# Date: November, 2022                                                        #
###############################################################################

# DESCRIPTION: This script is the overall executable bash script that executes
# all other bash and python scripts in the process, tracks the time and writes
# all progress to the logfile. Parallelization is done by invocing the scripts
# several times on subsets of the unzipped raw data.

# Dependencies:
# - `agrep` availabe via e.g. https://github.com/Wikinaut/agrep
# - `tre-agrep` available via e.g. https://wiki.ubuntuusers.de/tre-agrep/
# - `python`

# import config
CONFIG=${1}
source "${CONFIG}"

# export config and config variables to environment
export CONFIG=${CONFIG}
export ROOT_DIR=${ROOT_DIR}
export TMP_DIR=${TMP_DIR}
export OUT_DIR=${OUT_DIR}
export RES_DIR=${RES_DIR}
export PLT_DIR=${PLT_DIR}

# check folders
bash ${SCRIPT_0}

# create a logfile with date and time and export it
NOW=$(date +"%F_%H-%M")
export LOGFILE="${LOG_DIR}/log-$NOW.log"
touch "${LOGFILE}"

# create read file that stores number of reads
rm -r "${READ_FILE}"
printf "READ_FILE of ${PROJ_NAME}\n" > "${READ_FILE}"

# write header with parameters to log file
printf "%s\n" \
  "########################" \
  "Logfile of ${PROJ_NAME}" \
  "########################" \
  "Date and time: $(timestamp)" \
  "Output path: ${OUT_DIR}" \
  "Result path: ${RES_DIR}" \
  "Name of logfile: ${LOGFILE}" \
  "########################" \
  "" >> "${LOGFILE}"

# time of overall script start
START_0=$(date +%s)


################################## script 1 ##################################
# print current script to STDOUT and logfile (same for all other scripts)
printf "$(timestamp): starting script ${SCRIPT_1} ...\n" | tee -a "${LOGFILE}"


# start time (same for all other scripts)
START=$(date +%s)

# execute script (same for all other scripts)
bash ${SCRIPT_1}

# end time (same for all other scripts)
END=$(date +%s)

# print runtime to STDOUT and logfile
printf "$(timestamp): script ${SCRIPT_1} successfully executed\n" | tee -a "${LOGFILE}"
printf $((END-START)) | awk '{printf "runtime: %02d:%02d:%02d\n\n", int($1/3600), int($1/60)-int($1/3600)*60, int($1%60)}' | tee -a "${LOGFILE}"


################################## script 2 ##################################
printf "$(timestamp): starting script ${SCRIPT_2} ...\n" | tee -a "${LOGFILE}"
START=$(date +%s)

# execute script
Rscript ${SCRIPT_2}

END=$(date +%s)
printf "$(timestamp): script ${SCRIPT_2} successfully executed\n" | tee -a "${LOGFILE}"
printf $((END-START)) | awk '{printf "runtime: %02d:%02d:%02d\n\n", int($1/3600), int($1/60)-int($1/3600)*60, int($1%60)}' | tee -a "${LOGFILE}"


################################## script 3 ##################################
printf "$(timestamp): starting script ${SCRIPT_3} ...\n" | tee -a "${LOGFILE}"
START=$(date +%s)

# execute script
bash ${SCRIPT_3}

END=$(date +%s)
printf "$(timestamp): script ${SCRIPT_3} successfully executed\n" | tee -a "${LOGFILE}"
printf $((END-START)) | awk '{printf "runtime: %02d:%02d:%02d\n\n", int($1/3600), int($1/60)-int($1/3600)*60, int($1%60)}' | tee -a "${LOGFILE}"


################################## script 4 ##################################
printf "$(timestamp): starting script ${SCRIPT_4} ...\n" | tee -a "${LOGFILE}"
START=$(date +%s)

# execute script
bash ${SCRIPT_4}

END=$(date +%s)
printf "$(timestamp): script ${SCRIPT_4} successfully executed\n" | tee -a "${LOGFILE}"
printf $((END-START)) | awk '{printf "runtime: %02d:%02d:%02d\n\n", int($1/3600), int($1/60)-int($1/3600)*60, int($1%60)}' | tee -a "${LOGFILE}"


################################## script 5 ##################################
printf "$(timestamp): starting script ${SCRIPT_5} ...\n" | tee -a "${LOGFILE}"
START=$(date +%s)

# execute script
Rscript ${SCRIPT_5}

END=$(date +%s)
printf "$(timestamp): script ${SCRIPT_5} successfully executed\n" | tee -a "${LOGFILE}"
printf $((END-START)) | awk '{printf "runtime: %02d:%02d:%02d\n\n", int($1/3600), int($1/60)-int($1/3600)*60, int($1%60)}' | tee -a "${LOGFILE}"


################################## script 6 ##################################
printf "$(timestamp): starting script ${SCRIPT_6} ...\n" | tee -a "${LOGFILE}"
START=$(date +%s)

# execute script
bash ${SCRIPT_6}

END=$(date +%s)
printf "$(timestamp): script ${SCRIPT_6} successfully executed\n" | tee -a "${LOGFILE}"
printf $((END-START)) | awk '{printf "runtime: %02d:%02d:%02d\n\n", int($1/3600), int($1/60)-int($1/3600)*60, int($1%60)}' | tee -a "${LOGFILE}"


################################## script 7 ##################################
printf "$(timestamp): starting script ${SCRIPT_7} ...\n" | tee -a "${LOGFILE}"
START=$(date +%s)

# execute script
Rscript ${SCRIPT_7}

END=$(date +%s)
printf "$(timestamp): script ${SCRIPT_7} successfully executed\n" | tee -a "${LOGFILE}"
printf $((END-START)) | awk '{printf "runtime: %02d:%02d:%02d\n\n", int($1/3600), int($1/60)-int($1/3600)*60, int($1%60)}' | tee -a "${LOGFILE}"


############################# final time function #############################
# time of overall script end
END_0=$(date +%s)

# print total runtime
printf "########################\n" | tee -a "${LOGFILE}"
printf "Script successfully executed!\n" | tee -a "${LOGFILE}"
printf $((END_0-START_0)) | awk '{printf "Total runtime: %02d:%02d:%02d\n", int($1/3600), int($1/60)-int($1/3600)*60, int($1%60)}' | tee -a "${LOGFILE}"
printf "########################\n" | tee -a "${LOGFILE}"

# print end statement
printf "The End!\n"
