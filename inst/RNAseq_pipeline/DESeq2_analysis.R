#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Test if there are 12 args: if not, return an error
if (length(args)!=8) {
  stop("Please specify:\n
       [required] 1/ A comma-separated list of count (ref genome)\n
       [required] 2/ A comma-separated list of sample names \n
       [required] 3/ A comma-separated list of conditions \n
       [required] 4/ A comma-separated list of controls \n
       [required] 5/ dds output folder \n
       [required] 6/ FC tables output folder \n
       [required] 7/ PDF output folder \n
       [required] 8/ Experiment \n")
}

# Load libraries
suppressMessages(library(vlfunctions, warn.conflicts = FALSE))
suppressMessages(library(DESeq2, warn.conflicts = FALSE))

# Tests ----
# counts <- c("db/count_tables/RNAseq/Pprc1.depl/Pprc1.AID_IAA_0h_Pprc1.depl_rep1_mm10_counts.txt",
#             "db/count_tables/RNAseq/Pprc1.depl/Pprc1.AID_IAA_0h_Pprc1.depl_rep2_mm10_counts.txt",
#             "db/count_tables/RNAseq/Pprc1.depl/Pprc1.AID_IAA_1h_Pprc1.depl_rep1_mm10_counts.txt",
#             "db/count_tables/RNAseq/Pprc1.depl/Pprc1.AID_IAA_1h_Pprc1.depl_rep2_mm10_counts.txt",
#             "db/count_tables/RNAseq/Pprc1.depl/Pprc1.AID_IAA_3h_Pprc1.depl_rep1_mm10_counts.txt",
#             "db/count_tables/RNAseq/Pprc1.depl/Pprc1.AID_IAA_3h_Pprc1.depl_rep2_mm10_counts.txt",
#             "db/count_tables/RNAseq/Pprc1.depl/Pprc1.AID_IAA_6h_Pprc1.depl_rep1_mm10_counts.txt",
#             "db/count_tables/RNAseq/Pprc1.depl/Pprc1.AID_IAA_6h_Pprc1.depl_rep2_mm10_counts.txt")
# names <- c("control_rep1","control_rep2","IAA_1h_rep1","IAA_1h_rep2","IAA_3h_rep1","IAA_3h_rep2","IAA_6h_rep1","IAA_6h_rep2")
# conditions <- c("control","control","IAA_1h","IAA_1h","IAA_3h","IAA_3h","IAA_6h","IAA_6h")
# controls <- c("control","control","control","control","control","control","control","control")
# dds_output_folder <- "db/dds/RNAseq/"
# FC_output_folder <- "db/FC_tables/RNAseq/"
# PDF_output_folder <- "pdf/RNAseq/"
# experiment <- "Pprc1.depl"

# Parse variables ----
counts <- unlist(tstrsplit(args[1], ","))
names <- unlist(tstrsplit(args[2], ","))
names <- gsub("-", ".", names) # Names and conditions do not tolerate "-"
conditions <- unlist(tstrsplit(args[3], ","))
conditions <- gsub("-", ".", conditions) # Names and conditions do not tolerate "-"
controls <- unlist(tstrsplit(args[4], ","))
controls <- gsub("-", ".", controls) # Names and conditions do not tolerate "-"
dds_output_folder <- args[5]
FC_output_folder <- args[6]
PDF_output_folder <- args[7]
experiment <- args[8]

# Import data ----
dat <- lapply(counts, fread)
names(dat) <- names
dat <- rbindlist(dat, idcol = "condition")
dat[, condition:= factor(condition, unique(condition))]
DF <- dcast(dat, gene_id~condition, value.var = "count")
DF <- data.frame(DF[, -1], row.names = DF$gene_id)

# Remove low count reads ----
DF <- DF[rowSums(DF >= 10) >= 2,]

# SampleTable ----
sampleTable <- data.frame(condition = conditions,
                          row.names = names)

# DESeq2 analysis ----
print(paste("Start DESeq2 analysis normalization"))
# Create DESeq2 object ----
dds <- DESeqDataSetFromMatrix(countData = DF,
                              colData = sampleTable,
                              design = ~ condition)

# Compute model and save object ----
dds <- DESeq(dds)
saveRDS(dds,
        paste0(dds_output_folder, experiment, "/", experiment, "_DESeq2.dds"))

# Import combinations ----
cmb <- data.table(V1= conditions,
                  V2= controls)
cmb <- unique(cmb[V1!=V2])

# Compute FC and plot MA plots ----
cmb[, {
  res <- results(dds,
                 contrast = c("condition", V1, V2))
  res <- as.data.frame(res)
  res <- as.data.table(res, keep.rownames = T)
  res[, diff:= fcase(padj<0.05 & log2FoldChange>log2(1.5), "Up-regulated",
                     padj<0.05 & log2FoldChange<(-log2(1.5)), "Down-regulated",
                     default = "Unaffected")]
  # FC file ----
  outputFile <- paste0(FC_output_folder, experiment, "/", experiment, "_", V1, "_vs_", V2, "_DESeq2.txt")
  fwrite(res,
         outputFile,
         sep="\t",
         na = NA)
  print(paste(outputFile, "saved"))
  sumup <- table(res$diff)
  print(paste(names(sumup), "= ", sumup, collapse = " | "))
  
  # MA plot ----
  outputPdf <- paste0(PDF_output_folder, experiment, "/MA_plots/", experiment, "_", V1, "_vs_", V2, "_DESeq2_MA_plot.pdf")
  pdf(outputPdf, 4, 3)
  par(mai= c(.9,1.5,.9,1.3),
      cex.axis= 6/12,
      cex.lab= 7/12,
      las= 1,
      tcl= -0.1,
      mgp= c(.8, 0.25, 0),
      font.main= 1,
      cex.main= 9/12)
  vl_MAplot(DESeq2_dat= res,
            padj.cutoff = 0.05,
            log2FC.cutoff = log2(1.5),
            main= paste(V1, "vs.", V2))
  dev.off()
  print(paste(outputPdf, "saved"))
}, .(V1, V2)]

