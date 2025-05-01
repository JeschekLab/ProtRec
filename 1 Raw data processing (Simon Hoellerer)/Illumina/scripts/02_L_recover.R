# import variables
TMP_DIR <- Sys.getenv("TMP_DIR")

# change directory
setwd(TMP_DIR)

# load libraries
library(matrixStats)

# find best match to barcode and discriminator versions
for (FC in 1:2) {
  print(paste0("Reading data of flow cell ", FC, " ..."))
  # read distance as data frame
  df <- read.table(paste0("01_data_distance_all_", FC, ".tmp"), sep = "\t", header = F, colClasses = "integer")

  # rename according to their state and time point
  names(df) <- c(paste0("N", 1:6), paste0("F", 1:6))

  # find the minimum distance per row
  print(paste0("Finding minimum of flow cell ", FC, " ..."))
  df$min <- rowMins(as.matrix(df[, 1:12]))

  # find which column contains the minimum
  print(paste0("Finding which minimum of flow cell ", FC, " ..."))
  df$which_min <- max.col(-df[, 1:12])

  # count the appearance of the minimum
  print(paste0("Counting minimum of flow cell ", FC, " ..."))
  df$n_min <- rowSums(df[, 1:12] == df$min)

  # make printable column
  print_names <- names(df)[df$which_min]

  # ignore lines that have their minimum more than once
  print_names[df$n_min != 1] <- 0

  # ignore lines that have a minimum distance greater than 3
  print_names[df$min > 3] <- 0

  # make writable file
  write.table(print_names, file = paste0("02_nearest_neighbour_", FC, ".txt"), sep = "\t", quote = F,
    col.names = F, row.names = F)

  print(paste0("Done with FC", FC, "!"))
}
