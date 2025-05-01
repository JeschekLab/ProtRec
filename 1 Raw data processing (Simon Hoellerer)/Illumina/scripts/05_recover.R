# import variables
TMP_DIR <- Sys.getenv("TMP_DIR")

# change directory
setwd(TMP_DIR)

# load libraries
library(matrixStats)

# read data
for (state in c("N", "F")) {
  for (i in 1:6) {
    print(paste0("current state: ", state, " L", i))
    # read distance as data frame
    df <- read.table(paste0("04_L", i, "_", state, "_distance_all.tmp"), sep = "\t", header = F, colClasses = "integer")
     
    # rename according to their state and time point
    names(df) <- paste0("R", 1:6)

    # find the minimum distance per row
    df$min <- rowMins(as.matrix(df[, 1:6]))

    # find which column contains the minimum
    df$which_min <- max.col(-df[, 1:6])

    # count the appearance of the minimum
    df$n_min <- rowSums(df[, 1:6] == df$min)

    # make printable column
    print_names <- names(df)[df$which_min]

    # ignore lines that have their minimum more than once
    print_names[df$n_min != 1] <- 0

    # ignore lines that have a minimum distance greater than 3
    print_names[df$min > 3] <- 0

    # make writable file
    write.table(print_names, file = paste0("05_nearest_neighbour_L", i, "_", state, ".txt"), sep = "\t", quote = F,
      col.names = F, row.names = F)
  }
}
