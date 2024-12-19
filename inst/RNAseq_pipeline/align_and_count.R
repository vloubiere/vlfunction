#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least 2 args: if not, return an error
if (length(args)!=8) {
  stop("Please specify:\n
       [required] 1/ A comma-separated list of fq1 files \n
       [required] 2/ A comma-separated list of fq2 files \n
       [required] 3/ Is the data paired-end?
       [required] 4/ subreadr index prefix \n
       [required] 5/ Output bam path \n
       [required] 6/ Path to the gtf file that was used to generate the subreadr index \n
       [required] 7/ Output statistics file \n
       [required] 8/ Output count file \n")
}

suppressMessages(library(data.table))
suppressMessages(library(Rsubread))

# Tests ----
# fq1 <- "/scratch/stark/vloubiere/RNAseq/fq/223TCHLT1_2_R18054_fqfull_20241108_TCATTC_ACGTCCTG_1.fq.gz"
# fq2 <- "/scratch/stark/vloubiere/RNAseq/fq/223TCHLT1_2_R18054_fqfull_20241108_TCATTC_ACGTCCTG_2.fq.gz"
# paired <- TRUE
# idx <- "/groups/stark/vloubiere/genomes/Mus_musculus/subreadr_mm10/subreadr_mm10_index"
# bam <- "/scratch/stark/vloubiere/RNAseq/bam/test.bam"
# gtf <- "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/gencode.vM25.basic.annotation.gtf.gz"
# stats_file <- "db/counts/RNAseq/test_counts.txt"
# count_file <- "db/count_tables/RNAseq/test_counts.txt"

# Parse arguments ----
fq1 <- unlist(tstrsplit(args[1], ","))
fq2 <- unlist(tstrsplit(args[2], ","))
paired <- as.logical(args[3])
idx <- args[4]
bam <- args[5]
gtf <- args[6]
stats_file <- args[7]
count_file <- args[8]

# Catenate fq files if necessary
if(length(fq1)>1)
{
  temp_file <- tempfile(fileext = ".txt")
  system(paste("cat", paste(fq1, collapse = " "), ">", temp_file))
  fq1 <- temp_file
}
if(paired && length(fq2)>1)
{
  temp_file <- tempfile(fileext = ".txt")
  system(paste("cat", paste(fq2, collapse = " "), ">", temp_file))
  fq2 <- temp_file
}

# Align ----
Rsubread::align(index= idx,
                readfile1= fq1,
                readfile2= ifelse(paired, fq2, NULL),
                type= "rna",
                input_format = "gzFASTQ",
                output_format = "BAM",
                maxMismatches = 3,
                nthreads = data.table::getDTthreads()-1,
                unique = T,
                output_file= bam)

# Count ----
.c <- Rsubread::featureCounts(bam,
                              annot.ext= gtf,
                              isGTFAnnotationFile = TRUE,
                              GTF.featureType= "exon",
                              GTF.attrType= "gene_id",
                              GTF.attrType.extra= "gene_name",
                              isPairedEnd = paired,
                              nthreads = data.table::getDTthreads()-1)

# Save statistics ----
aln_stats <- fread(paste0(bam, ".summary"))
setnames(aln_stats, c("feature", "count"))
assign_stats <- .c$stat
setnames(assign_stats, c("feature", "count"))
stats <- rbind(aln_stats, assign_stats)
fwrite(stats,
       stats_file,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote= F)

# Save count table ----
counts <- data.table(gene_id= as.data.table(.c$annotation)[, paste0(GeneID, "__", gene_name, "__", Length)],
                     count= .c[[1]][,1])
fwrite(counts,
       count_file,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote= F)
