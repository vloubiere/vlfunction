#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least 2 args: if not, return an error
if (length(args)!=2) {
  stop("Please specify:\n
       [required] 1/ Input bam file \n
       [required] 2/ Output bw file \n")
}

suppressMessages(library(vlfunctions, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomicRanges, warn.conflicts = FALSE))

# # Test ----
# bam_file <- "/scratch/stark/vloubiere/RNAseq/bam/test.bam"
# bw_file <- "/groups/stark/vloubiere/projects/vl_pipelines/db/bw/RNAseq/test.bw"

# Parse args ----
bam_file <- args[1]
bw_file <- args[2]

# Variables test ----
bed <- vl_importBam(file = bam_file,
                    sel = c("qname", "flag", "rname", "pos", "qwidth"),
                    col.names = c("readID", "samFlag", "seqnames", "start", "width"))

# Identify unambiguously mapped paired-end reads ----
bed[, check:= (bitwAnd(samFlag, 0x4) == 0) & (bitwAnd(samFlag, 0x2) != 0)]
bed <- bed[(check)]

# Extend reads ----
bed[, end:= start+width-1]
setorderv(bed,
          c("seqnames", "start", "end"))

# Compute coverage and export ----
gr <- GenomicRanges::GRanges(bed)
cov <- GenomicRanges::coverage(gr)
rtracklayer::export(cov,
                    bw_file)
