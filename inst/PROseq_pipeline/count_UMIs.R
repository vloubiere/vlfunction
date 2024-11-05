#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least 2 args: if not, return an error
if (length(args)!=2) {
  stop("Please specify:\n
       [required] 1/ Aligned bam file containing UMIs\n
       [required] 2/ Output count file (.txt)\n")
}

suppressMessages(library(vlfunctions, warn.conflicts = FALSE))
suppressMessages(library(Rsamtools, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomicRanges, warn.conflicts = FALSE))
suppressMessages(library(stringdist, warn.conflicts = FALSE))

# Variables ----
# bam <- "/scratch/stark/vloubiere/PROSeq_pipeline/bam/truncated.bam"
# output <- "db/umi_counts/test.txt"
bam <- args[1]
output <- args[2]

# Import data ----
dat <- vl_importBam(bam,
                    sel= c("qname", "rname", "pos", "strand", "qwidth", "mapq"),
                    col.names = c("read", "seqnames", "start", "strand", "width", "mapq"))
dat[, strand:= as.character(strand)]
dat[!is.na(strand), start:= ifelse(strand=="+", start, start+width-1)]
dat[!is.na(strand), strand:= ifelse(strand=="+", "-", "+")] # Invert strand
dat[, coor:= paste0(seqnames, ":", start, ":", strand)]

# Extract UMIs ----
dat[!is.na(seqnames), UMI:= gsub(".*:", "", read)]
dat[is.na(seqnames), UMI:= "GGGGGGGGGG"] # Unaligned reads (not collapsed, keep for total reads)
dat <- dat[, .(umi_N= .N), .(coor, UMI)]
dat[, total_counts:= sum(umi_N), coor]
setorderv(dat, "umi_N", order = -1)

# Check whether UMI might be collapsed ----
dat[, collapsed:= T, coor]
dat[, idx:= .I]
for(i in 1:10)
{
  dat[, check:= idx[1], .(coor, gsub(paste0("^(.{", i-1, "})."), "\\1", UMI))]
  potentialDup <- unique(dat[(check<idx), c(check, idx)])
  dat[potentialDup, collapsed:= FALSE]
  print(i)
}
dat$idx <- NULL
paste0(sum(dat$collapsed), " / ", nrow(dat), " pre-collapsed")

# UMI collapsing (>1 diff) ----
while(any(!dat$collapsed))
{
  dat[!(collapsed), c("collapsed", "UMI"):= {
    coll <- stringdist(UMI[1],
                       UMI,
                       method="hamming",
                       nthread= getDTthreads()-1)<=1
    UMI[coll] <- UMI[1]
    .(coll, UMI)
  }, coor]
}

# Final collapsing ----
dat <- unique(dat[, .(coor, total_counts, UMI)])
dat <- dat[, .(umi_counts= .N), .(coor, total_counts)]
dat[coor=="NA:NA:NA", umi_counts:= NA]

# Save output ----
fwrite(dat,
       output,
       sep= "\t",
       quote= F,
       na= NA)
