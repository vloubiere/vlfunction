#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# test if there is at least 2 args: if not, return an error
if (length(args)!=2) {
  stop("Please specify:\n
       [required] 1/ Count file \n
       [required] 2/ Output stat file")
}
require(data.table)

# Compute statistics and spikein sizeFactor ----
dat <- fread(args[1])
stats <- data.table(total= sum(dat$total_counts),
                    mapped= sum(dat[coor!="NA:NA:NA", total_counts]), # Remove unmapped
                    umi_counts= sum(dat[coor!="NA:NA:NA", umi_counts]))
fwrite(stats,
       args[2],
       sep= "\t",
       na= NA)