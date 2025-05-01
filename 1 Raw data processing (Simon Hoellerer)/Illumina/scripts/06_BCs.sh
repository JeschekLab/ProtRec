#!/bin/bash

# import config
source "${CONFIG}"

# change directory
cd "${TMP_DIR}"

# combine data and aligned read
for state in "N" "F"; do
  for t in {1..6}; do
    paste \
      "03_${state}_L${t}.txt" \
      "05_nearest_neighbour_L${t}_${state}.txt" \
      > "06_idx_L${t}_${state}.tmp" &
  done
done
wait

# split into barcode combinations
for state in "N" "F"; do
  for L in {1..6}; do
    for R in {1..6}; do
      awk -v pat=${R} '$2 ~ pat' "06_idx_L${L}_${state}.tmp" \
        > "06_L${L}_R${R}_${state}.txt" &
    done
  done
done
wait

# trim indices off reads
for L in {1..6}; do
  for state in "N" "F"; do
    for R in {1..6}; do
      let Rstart=(R+9)
      printf "Trim barcodes for 06_L${L}_R${R}_${state} at ${Rstart} ...\n"
      cat "06_L${L}_R${R}_${state}.txt" \
        | awk '{print substr($1, '${Rstart}', 57)}' \
        > "06_L${L}_R${R}_${state}_trimmed.txt" &
    done
  done
  wait
done

# sort files
for state in "N" "F"; do
  for L in {1..6}; do
    for R in {1..6}; do
      printf "Sorting 06_L${L}_R${R}_${state}_trimmed.txt ...\n"
      sort --parallel=12 -S 80G "06_L${L}_R${R}_${state}_trimmed.txt" \
        > "06_L${L}_R${R}_${state}_sorted.tmp" &
    done
    wait
  done
done

# count lines
printf "\nSorted and trimmed\n" >> "${READ_FILE}"
for state in "N" "F"; do
  for L in {1..6}; do
    for R in {1..6}; do
      wc -l "06_L${L}_R${R}_${state}_sorted.tmp" \
        | tr '/' '\t' | awk 'BEGIN {OFS = "\t"} {print $2, $1}' \
        >> "${READ_FILE}"
      printf "Done with 06_L${L}_R${R}_${state}_sorted.tmp ...\n"
    done
  done
done

# count uniques
for state in "N" "F"; do
  for L in {1..6}; do
    for R in {1..6}; do
      printf "Uniques in 07_L${L}_R${R}_${state}_sorted.tmp ...\n"
      uniq -c "06_L${L}_R${R}_${state}_sorted.tmp" \
        | awk 'BEGIN {OFS = "\t"} {print $2, $1}' \
        > "06_L${L}_R${R}_${state}_uniq.tmp" &
      sleep 0.5
    done
  done
done
wait

# move to output directory
for state in "N" "F"; do
  for L in {1..6}; do
    for R in {1..6}; do
      printf "Moving 06_L${L}_R${R}_${state}_uniq.tmp ...\n"
      cp "06_L${L}_R${R}_${state}_uniq.tmp" "${OUT_DIR}/L${L}_R${R}_${state}.txt"
      sleep 0.5
    done
  done
done

# done!
