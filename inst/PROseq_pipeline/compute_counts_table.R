#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least 2 args: if not, return an error
if (length(args)!=3) {
  stop("Please specify:\n
       [required] 1/ Reference genome umi count file \n
       [required] 2/ A comma-separated list of annotation files \n
       [required] 3/ A comma-separated list of output file names \n")
}

require(data.table)

# Tests ----
# umi_count <- "db/counts/PROseq/HCFC1/AID-Hcfc1-cl4_0hrIAA_HCFC1_rep1_mm10_counts.txt"
# annotations <- c("/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_promoters.rds",
#                  "/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_genebody.rds",
#                  "/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_transcript.rds")
# outputFiles <- c("db/count_tables/PROseq/HCFC1/promoter/AID-Hcfc1-cl4_0hrIAA_HCFC1_rep1_mm10_promoter_counts.txt",
#                 "db/count_tables/PROseq/HCFC1/gene_body/AID-Hcfc1-cl4_0hrIAA_HCFC1_rep1_mm10_gene_body_counts.txt",
#                 "db/count_tables/PROseq/HCFC1/transcript/AID-Hcfc1-cl4_0hrIAA_HCFC1_rep1_mm10_transcript_counts.txt")

# Args ----
umi_count <- args[1]
annotations <- unlist(tstrsplit(args[2], ","))
outputFiles <- unlist(tstrsplit(args[3], ","))

# Import annotation ----
annot <- data.table(file= annotations,
                    output= outputFiles)
annot <- annot[, readRDS(file), output]

# Import UMI ----
dat <- fread(umi_count)
dat[, c("seqnames", "start", "strand"):= tstrsplit(coor, ":", type.convert = T)]

# Compute counts ----
annot$count <- dat[annot, sum(umi_counts, na.rm= T), on= c("seqnames", "start>=start", "start<=end", "strand"), .EACHI]$V1
annot[, {
  fwrite(.SD[, .(ID= cluster.id, count)],
         output,
         sep= "\t",
         na= NA)
}, output]
