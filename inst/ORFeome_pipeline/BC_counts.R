#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

############----------------------------------------------------##############
########################       Counts per BC        ##########################
############----------------------------------------------------##############

# Test if arguments provided
if (length(args)!=4) {
  stop("Please specify:\n
       [required] 1/ Aligned BCs bam\n
       [required] 2/ Dictionary file (.rds)\n
       [required] 3/ Output stats file (.txt)\n
       [required] 4/ Output counts file (.txt)\n")
}

require(Rsamtools)
require(rtracklayer)
require(GenomicRanges)
require(data.table)

# Test ----
# bam <- "/scratch/stark/vloubiere/ORFeome/bam/apoptosisFasL_FasL_NA_A549_rep1_lib200.bam"
# BC <- "/groups/stark/vloubiere/projects/viralORF_tomas/db/dictionary/lib200_merged_dictionary.rds"
# output_stats <- "db/alignment_stats/ORFeome/apoptosisFasL_FasL_NA_A549_rep1_lib200_stats.txt"
# output_counts <- "db/counts/ORFeome/apoptosisFasL_FasL_NA_A549_rep1_lib200_counts.txt"

# Variables ----
bam <- args[1]
BC <- args[2]
output_stats <- args[3]
output_counts <- args[4]

# Import bam ----
.c <- Rsamtools::scanBam(bam,
                         param = Rsamtools::ScanBamParam(what= c("qname", "rname", "mapq", "strand", "seq")))
.c <- .c[[1]]
.c$seq <- as.character(.c$seq)
.c <- as.data.table(.c)
.c[, rname:= as.character(rname)]
.c[, strand:= as.character(strand)]

# Alignment statistics ----
stats <- .c[, .(total= .N,
                aligned= sum(!is.na(rname)),
                "mapq>=30"= sum(mapq>=30, na.rm = T))]
fwrite(stats,
       output_stats,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)

# Select aligned reads with mapq>30 ----
.c <- .c[!is.na(rname) & mapq>=30]
.c <- .c[, .(count= sum(mapq>=30, na.rm = T)), rname]

# Import dictionary ----
dic <- readRDS(BC)
setnames(dic, "bcID", "ID")
dic <- unique(dic[, .(ID, BC)])
dic[.c, count:= i.count, on= "ID==rname"]
dic[is.na(count), count:= 0]

# Save ----
fwrite(dic,
       output_counts,
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)