#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least 2 args: if not, return an error
if (length(args)!=3) {
  stop("Please specify:\n
       [required] 1/ Reference genome umi count file \n
       [required] 2/ A .rds annotation file \n
       [required] 3/ Output file name \n")
}

suppressMessages(library(data.table))

# Tests ----
# umi_count <- "db/counts/PROseq/HCFC1/AID-Hcfc1-cl4_0hrIAA_HCFC1_rep1_mm10_counts.txt"
# annotation <- "/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_promoters.rds"
# outputFile <- "db/count_tables/PROseq/HCFC1/promoter/AID-Hcfc1-cl4_0hrIAA_HCFC1_rep1_mm10_promoter_counts.txt"

# Args ----
umi_count <- args[1]
annotation <- args[2]
outputFile <- args[3]

# Import annotation ----
annot <- readRDS(annotation)

# Import UMI ----
dat <- fread(umi_count)
dat[, c("seqnames", "start", "strand"):= tstrsplit(coor, ":", type.convert = T)]

# Compute counts ----
annot$count <- dat[annot, sum(umi_counts, na.rm= T), on= c("seqnames", "start>=start", "start<=end", "strand"), .EACHI]$V1
fwrite(annot[, .(ID= cluster.id, count)],
       outputFile,
       sep= "\t",
       na= NA)