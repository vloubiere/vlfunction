#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least 2 args: if not, return an error
if (length(args)!=2) {
  stop("Please specify:\n
       [required] 1/ UMI counts file\n
       [required] 2/ Output prefix (.ps.bw; .ns.bw)\n")
}

suppressMessages(library(vlfunctions, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomicRanges, warn.conflicts = FALSE))

# Variables test
# countsFile <- "db/umi_counts/HCFC1/AID-Hcfc1-cl17_0hrIAA_rep1_mm10_counts.txt"

# Variables ----
countsFile <- args[1]
outputPrefix <- args[2]
  
# Import umi counts ----
counts <- fread(countsFile)
counts <- counts[coor!="NA:NA:NA"]
counts[, c("seqnames", "start", "strand"):= tstrsplit(coor, ":", type.convert = T)]

# Expand counts ----
counts <- counts[rep(seq(.N), umi_counts), .(seqnames, start, strand)]
counts[, end:= start]
counts <- vl_resizeBed(counts, "start", 0, 9) # Resize to improve visualization
setorderv(counts, c("seqnames", "start", "end"))

# Compute coverage and export ----
counts[, {
  gr <- GenomicRanges::GRanges(.SD)
  cov <- GenomicRanges::coverage(gr)
  outputFile <- paste0(outputPrefix, ifelse(strand=="+", ".ps", ".ns"), ".bw")
  rtracklayer::export(cov,
                      outputFile)
  print(outputFile)
}, strand]
