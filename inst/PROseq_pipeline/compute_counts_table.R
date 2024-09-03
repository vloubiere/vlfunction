#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least 2 args: if not, return an error
if (length(args)!=6) {
  stop("Please specify:\n
       [required] 1/ Reference genome count file \n
       [required] 2/ Spike-in count file \n
       [required] 3/ Output folder \n
       [required] 4/ .rds file containing promoter annotations \n
       [required] 5/ .rds file containing gene body annotations \n
       [required] 6/ .rds file containing transcript annotations \n")
}

require(data.table)

# Tests ----
# umi_count <- "db/umi_counts/AID-Hcfc1-cl4_0hrIAA_rep1_mm10_counts.txt"
# spike <- "db/umi_counts/AID-Hcfc1-cl4_0hrIAA_rep1_dm3_spikein_counts.txt"
# outputFolder <- "db/count_tables/HCFC1/"
# annotations <- c("/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_promoters.rds",
#                  "/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_genebody.rds",
#                  "/groups/stark/vloubiere/projects/PROseq_pipeline/db/annotations/mm10_transcript.rds")

# Args ----
umi_count <- args[1]
spike <- args[2]
outputFolder <- args[3]
annotations <-  args[4:6]

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

# Compute stats ----
stats <- data.table(total= sum(dat$total_counts),
                    mapped= sum(dat[coor!="NA:NA:NA", total_counts]), # Remove unmapped
                    umi_counts= sum(dat[coor!="NA:NA:NA", umi_counts]))

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

# Compute statistics and spikein sizeFactor ----
datSpike <- fread(spike)
datSpike <- datSpike[coor!="NA:NA:NA"] # Remove unmapped
datSpike <- data.table(total_spikeIn= sum(datSpike$total_counts),
                       umi_counts_spikeIn= sum(datSpike$umi_counts, na.rm= T))
stats <- cbind(stats, datSpike)
fwrite(stats,
       paste0(outputFolder, "/stats/", gsub("counts.txt$", "spikeIn_statistics.txt", basename(umi_count))),
       sep= "\t",
       na= NA)