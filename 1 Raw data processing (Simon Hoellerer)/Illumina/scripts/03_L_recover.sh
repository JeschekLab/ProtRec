#!/bin/bash

# import config
source "${CONFIG}"

# change directory
cd "${TMP_DIR}"

# combine remaining data and nearest BC/disc neighbour
for FC in {1..2}; do
  paste "01_data_rest_${FC}.txt" "02_nearest_neighbour_${FC}.txt" \
    > "03_data_rest_idx_${FC}.tmp" &
done
wait

# extract recovered reads to new files
for FC in {1..2}; do
  for state in "N" "F"; do
    for t in {1..6}; do
      awk -v pat="${state}${t}" '$3 ~ pat' "03_data_rest_idx_${FC}.tmp" \
        | awk '{print $2}' \
        > "03_recovered_${state}${t}_${FC}.tmp" &
    done
  done
  wait
done

# copy files to prevent overwriting
for FC in {1..2}; do
  for state in "N" "F"; do
    for t in {1..6}; do
      printf "${FC}_${state}_${t}\n"
      awk '{print $2}' "01_${state}_L${t}_${FC}.tmp" \
        > "03_${state}_L${t}_recovered_${FC}.tmp" &
    done
  done
  wait
done

# write not aligned reads to rest file
for FC in {1..2}; do
  awk '$3 ~ /0/' "03_data_rest_idx_${FC}.tmp" \
    | awk '{print $1, $2}' \
    > "03_data_unrecovered_${FC}.tmp" &
done
wait

# combine perfect and recovered reads
for FC in {1..2}; do
  for state in "N" "F"; do
    for t in {1..6}; do
      printf "${FC}_${state}_${t}\n"
      cat "03_${state}_L${t}_recovered_${FC}.tmp" \
        "03_recovered_${state}${t}_${FC}.tmp" \
        > "03_${state}_L${t}_${FC}.txt" &
    done
  done
done
wait

# combine both FCs
for FC in {1..2}; do
  for state in "N" "F"; do
    for t in {1..6}; do
      printf "${FC}_${state}_${t}\n"
      cat "03_${state}_L${t}_${FC}.txt" \
        >> "03_${state}_L${t}.txt" &
    done
  done
  wait
done

# count lines
printf "\nMapped left BC and discriminator\n" >> "${READ_FILE}"
for state in "N" "F"; do
  for t in {1..6}; do
    printf "${state}_${t}\n"
    wc -l "03_${state}_L${t}.txt" | tr '/' '\t' \
      | awk 'BEGIN {OFS = "\t"} {print $2, $1}' >> "${READ_FILE}"
  done
done

for FC in {1..2}; do
  wc -l "03_data_unrecovered_${FC}.tmp" | tr '/' '\t' \
    | awk 'BEGIN {OFS = "\t"} {print $2, $1}' >> "${READ_FILE}"
done
