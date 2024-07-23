#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

############ Assign insertions to closest downstream exon start ##############
# Takes as input a collapsed bam file, a path to gtf file containing non-first exons, a name for the output bed file and for the assignment table 
# 1/ Import bam file, compute insertion coordinates and collapse
# 2/ Save insertion as bed file (using provided name)
# 3/ Compute closest downstream non-first exon junction, both on the same and reverse strands
#    -> saved as: assignment_table"_same_strand.txt" & assignment_table"_rev_strand.txt"

# test if there is at least 2 args: if not, return an error
if (length(args)!=4) {
  stop("Please specify:\n
       [required] 1/ Collapsed bam file\n
       [required] 2/ Path to a gtf file containing non-first exons used for assignment\n
       [required] 3/ Output bed file (.bed)\n
       [required] 4/ Output assignment_table basename\n")
}

library(Rsamtools, warn.conflicts = FALSE)
library(rtracklayer, warn.conflicts = FALSE)
library(GenomicRanges, warn.conflicts = FALSE)
library(data.table, warn.conflicts = FALSE)

# For tests
# bam <- "db/bam_unique/Repressor1_input_rep1_q30_unique.bam"
# exons <- "db/gtf/exons_start_mm10.gtf"
# bed_file <- "test/test.bed"
# assignment_table <- "test/assignment_test"

# Variables ----
bam <- args[1]
exons <- args[2]
bed_file <- args[3]
assignment_table <- args[4]

# Import and parse bam file ----
.c <- Rsamtools::scanBam(bam,
                         param = Rsamtools::ScanBamParam(what=c("rname", "strand", "pos", "qwidth")))[[1]]
.c <- as.data.table(.c)
setnames(.c, c("seqnames", "strand", "start", "width"))
# Compute start
.c[strand=="+", start:= start-1] 
.c[strand=="-", start:= start+width]
.c[, strand:= ifelse(strand=="+", "-", "+")]# Reverse strand (iPCR)
.c[, end:= start]
# Remove non unique insertions
.c <- unique(.c)
.c <- na.omit(.c)
setorderv(.c, c("seqnames", "start", "end"))
# Save bed file
.c <- GenomicRanges::GRanges(.c)
export(.c, bed_file)
# Import non-first exons gtf file
exons <- rtracklayer::import(exons)
cols <- c("gene_id", "gene_name", "mgi_id", "exon_number", "exon_id")
cols <- cols[cols %in% names(mcols(exons))] # Human gtf does not contain mgi_id
mcols(exons) <- mcols(exons[, cols])
# For each strand, assign insertions to closest downstream, non-first exon ----
for(.strand in c("same_strand", "rev_strand"))
{
  curr <- .c
  if(.strand=="rev_strand")
    strand(curr) <- sapply(strand(curr), switch,  "-"="+", "+"="-", "*"="*")
  idx <- precede(curr, exons)
  mcols(curr) <- NULL
  mcols(curr) <- mcols(exons)[idx,]
  mcols(curr)$dist <- NA
  curr$dist[!is.na(idx)] <- distance(curr[!is.na(idx)], exons[idx[!is.na(idx)]])
  # SAVE
  curr <- as.data.table(curr)
  curr$width <- NULL
  fwrite(curr, 
         paste0(assignment_table, "_", .strand, ".txt"),
         col.names = T,
         quote= F,
         sep= "\t",
         row.names = F,
         na= NA)
}




