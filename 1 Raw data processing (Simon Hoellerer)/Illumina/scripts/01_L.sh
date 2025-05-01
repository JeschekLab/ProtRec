#!/bin/bash

# import config
source "${CONFIG}"

# change directory
cd "${TMP_DIR}"

# extract only first 34 bases of forward read (rest is irrelevant in this case)
paste \
  <(awk '{print substr($1, 1, 34)}' "${RAW_DIR}/data_all_LH.txt") \
  <(awk '{print $2}' "${RAW_DIR}/data_all_LH.txt") \
  > "00_data_cropped.txt"

# count lines
wc -l "00_data_cropped.txt" | tr '/' '\t' \
  | awk 'BEGIN {OFS = "\t"} {print $2, $1}' >> "${READ_FILE}"

# split file into 2 arbitrary flowcells and rename
split --lines=900000000 --suffix-length=1 --numeric-suffixes \
  "00_data_cropped.txt" 00_data_cropped_

mv "00_data_cropped_1" "00_data_cropped_2.txt"
mv "00_data_cropped_0" "00_data_cropped_1.txt"

# count lines
for FC in {1..2}; do
  wc -l "00_data_cropped_${FC}.txt" | tr '/' '\t' \
    | awk 'BEGIN {OFS = "\t"} {print $2, $1}' >> "${READ_FILE}"
done

# extract perfect L1 to L6 non-flipped and flipped into new files
for FC in {1..2}; do
  # loop through barcodes 1-6 with non-flipped discriminator
  for t in {1..6}; do
    let count=(t-1)
    pattern=${read_L_N[count]}
    printf "FC${FC}\tN${t}\tC${count}\t${pattern}\n"
    # find all perfectly matching reads
    awk -v pat=${pattern} '$1 ~ pat' "00_data_cropped_${FC}.txt" \
      > "01_N_L${t}_${FC}.tmp" &
  done
done
wait

for FC in {1..2}; do
  # loop through barcodes 1-6 with flipped discriminator
  for t in {1..6}; do
    let count=(t-1)
    pattern=${read_L_F[count]}
    printf "FC${FC}\tF${t}\tC${count}\t${pattern}\n"
    # find all perfectly matching reads
    awk -v pat=${pattern} '$1 ~ pat' "00_data_cropped_${FC}.txt" \
      > "01_F_L${t}_${FC}.tmp" &
  done
done
wait

# count lines
printf "\nCorrect left BC and discriminator\n" >> "${READ_FILE}"
for FC in {1..2}; do
  for state in "N" "F"; do
    for t in {1..6}; do
      printf "${state}_${t}\n"
      wc -l "01_${state}_L${t}_${FC}.tmp" | tr '/' '\t' \
        | awk 'BEGIN {OFS = "\t"} {print $2, $1}' >> "${READ_FILE}"
    done
  done
done

# write all reads that did not match any BC/disc version to a new file
for FC in {1..2}; do
  cat "00_data_cropped_${FC}.txt" \
    | awk -v pat="${read_L_N[0]}" '$1 !~ pat' \
    | awk -v pat="${read_L_N[1]}" '$1 !~ pat' \
    | awk -v pat="${read_L_N[2]}" '$1 !~ pat' \
    | awk -v pat="${read_L_N[3]}" '$1 !~ pat' \
    | awk -v pat="${read_L_N[4]}" '$1 !~ pat' \
    | awk -v pat="${read_L_N[5]}" '$1 !~ pat' \
    | awk -v pat="${read_L_F[0]}" '$1 !~ pat' \
    | awk -v pat="${read_L_F[1]}" '$1 !~ pat' \
    | awk -v pat="${read_L_F[2]}" '$1 !~ pat' \
    | awk -v pat="${read_L_F[3]}" '$1 !~ pat' \
    | awk -v pat="${read_L_F[4]}" '$1 !~ pat' \
    | awk -v pat="${read_L_F[5]}" '$1 !~ pat' \
    > "01_data_rest_${FC}.txt" &
  wait
done

# calculate distance between correct and wrong reads
for FC in {1..2}; do
    # loop through barcodes 1-6
  for t in {1..6}; do
    let count=(t-1)
    # map to non-flipped BC-discriminator combination
    cat "01_data_rest_${FC}.txt" | awk '{print $1}' \
      | tre --max-errors=34 --show-cost ${read_L_N[count]} \
      | tr ':' '\t' | awk '{print $1}' > "01_data_rest_distance_N${t}_${FC}.tmp" &
    # map to flipped BC-discriminator combination
    cat "01_data_rest_${FC}.txt" | awk '{print $1}' \
      | tre --max-errors=34 --show-cost ${read_L_F[count]} \
      | tr ':' '\t' | awk '{print $1}' > "01_data_rest_distance_F${t}_${FC}.tmp" &
  done
  wait
done

# paste together
for FC in {1..2}; do
  paste \
    "01_data_rest_distance_N1_${FC}.tmp" \
    "01_data_rest_distance_N2_${FC}.tmp" \
    "01_data_rest_distance_N3_${FC}.tmp" \
    "01_data_rest_distance_N4_${FC}.tmp" \
    "01_data_rest_distance_N5_${FC}.tmp" \
    "01_data_rest_distance_N6_${FC}.tmp" \
    "01_data_rest_distance_F1_${FC}.tmp" \
    "01_data_rest_distance_F2_${FC}.tmp" \
    "01_data_rest_distance_F3_${FC}.tmp" \
    "01_data_rest_distance_F4_${FC}.tmp" \
    "01_data_rest_distance_F5_${FC}.tmp" \
    "01_data_rest_distance_F6_${FC}.tmp" \
    > "01_data_distance_all_${FC}.tmp"
done

# count lines
for FC in {1..2}; do
  for t in {1..6}; do
    wc -l "01_data_rest_distance_N${t}_${FC}.tmp"| tr '/' '\t' \
        | awk 'BEGIN {OFS = "\t"} {print $2, $1}'
    wc -l "01_data_rest_distance_F${t}_${FC}.tmp"| tr '/' '\t' \
        | awk 'BEGIN {OFS = "\t"} {print $2, $1}'
  done
done
