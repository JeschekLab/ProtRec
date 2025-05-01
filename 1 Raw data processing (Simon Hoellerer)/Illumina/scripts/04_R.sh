#!/bin/bash

# import config
source "${CONFIG}"

# change directory
cd "${TMP_DIR}"

# write first 15 characters of reverse reads into new files
for state in "N" "F"; do
  for t in {1..6}; do
    awk '{print substr($1, 1, 15)}' "03_${state}_L${t}.txt" \
      > "04_substring_${state}_L${t}.tmp" &
  done
done
wait

# calculate distance between correct BC and substring of reads
for t in {1..6}; do
  printf "${t}\n"
  for state in "N" "F"; do
    for count in {0..5}; do
      let R=(count+1)
      tre --max-errors=15 --show-cost ${idx_R[count]} "04_substring_${state}_L${t}.tmp" \
      | tr ':' '\t' | awk '{print $1}' \
      > "04_distance_L${t}_${state}_R${R}.tmp" &
    done
  done
  wait
done

# paste into 12 files
for t in {1..6}; do
  for state in "N" "F"; do
    paste "04_distance_L${t}_${state}_R1.tmp" \
      "04_distance_L${t}_${state}_R2.tmp" \
      "04_distance_L${t}_${state}_R3.tmp" \
      "04_distance_L${t}_${state}_R4.tmp" \
      "04_distance_L${t}_${state}_R5.tmp" \
      "04_distance_L${t}_${state}_R6.tmp" \
      > "04_L${t}_${state}_distance_all.tmp" &
  done
done
wait
