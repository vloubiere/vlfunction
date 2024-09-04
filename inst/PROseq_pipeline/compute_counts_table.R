#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least 2 args: if not, return an error
if (length(args)!=5) {
  stop("Please specify:\n
       [required] 1/ Reference genome count file \n
       [required] 2/ Output folder \n
       [required] 3/ .rds file containing promoter annotations \n
       [required] 4/ .rds file containing gene body annotations \n
       [required] 5/ .rds file containing transcript annotations \n")
}

require(data.table)

# Tests ----
# umi_count <- "db/umi_counts/AID-Hcfc1-cl4_0hrIAA_rep1_mm10_counts.txt"
# outputFolder <- "db/count_tables/HCFC1/"
# annotations <- c("/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_promoters.rds",
#                  "/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_genebody.rds",
#                  "/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_transcript.rds")

# Args ----
umi_count <- args[1]
outputFolder <- args[2]
annotations <-  args[3:5]

# Create sub-directories ----
dir.create(paste0(outputFolder, "/promoter/"), showWarnings = F)
dir.create(paste0(outputFolder, "/gene_body/"), showWarnings = F)
dir.create(paste0(outputFolder, "/transcript/"), showWarnings = F)
dir.create(paste0(outputFolder, "/stats/"), showWarnings = F)

# Import annotations ----
annots <- data.table(feature= c("promoter", "gene_body", "transcript"),
                     file= annotations)
annots <- annots[, readRDS(file), feature]
annots[, seqnames:= as.character(seqnames)]
annots[, start:= as.integer(start)]
annots[, end:= as.integer(end)]
annots[, strand:= as.character(strand)]

# Import UMI ----
dat <- fread(umi_count)
dat[, c("seqnames", "start", "strand"):= tstrsplit(coor, ":", type.convert = T)]

# Compute counts ----
annots[, {
  .c <- .SD[, .(seqnames, start, end, strand)]
  .c <- data.table(ID= cluster.id,
                   count= dat[.c, sum(umi_counts, na.rm= T), on= c("seqnames", "start>=start", "start<=end", "strand"), .EACHI]$V1)
  outputFile <- paste0(outputFolder, "/", feature, "/", basename(umi_count))
  fwrite(.c,
         outputFile,
         sep= "\t",
         na= NA)
  print(outputFile)
}, feature]