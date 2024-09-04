#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(vlfunctions)
library(DESeq2)

# Test if there are 10 args: if not, return an error
if (length(args)!=10) {
  stop("Please specify:\n
       [required] 1/ A comma-separated list of count (ref genome)\n
       [required] 2/ A comma-separated list of read statistics (reference genome) \n
       [required] 3/ A comma-separated list of spike-in statistics \n
       [required] 4/ A comma-separated list of sample names \n
       [required] 5/ A comma-separated list of conditions \n
       [required] 6/ A comma-separated list of controls \n
       [required] 7/ .dds output prefix \n
       [required] 8/ dds statistics output \n
       [required] 9/ FC tables output folder \n
       [required] 10/ MAplot prefix \n")
}

# # Tests ----
# counts <- c("db/count_tables/PROseq/HCFC1/transcript/AID-Hcfc1-cl4_0hrIAA_HCFC1_rep1_mm10_counts.txt",
#             "db/count_tables/PROseq/HCFC1/transcript/AID-Hcfc1-cl4_3hrIAA_HCFC1_rep1_mm10_counts.txt",
#             "db/count_tables/PROseq/HCFC1/transcript/AID-Hcfc1-cl17_0hrIAA_HCFC1_rep1_mm10_counts.txt",
#             "db/count_tables/PROseq/HCFC1/transcript/AID-Hcfc1-cl17_3hrIAA_HCFC1_rep1_mm10_counts.txt")
# refStats <- c("db/count_tables/PROseq/HCFC1/stats/AID-Hcfc1-cl4_0hrIAA_HCFC1_rep1_mm10_statistics.txt",
#               "db/count_tables/PROseq/HCFC1/stats/AID-Hcfc1-cl4_3hrIAA_HCFC1_rep1_mm10_statistics.txt",
#               "db/count_tables/PROseq/HCFC1/stats/AID-Hcfc1-cl17_0hrIAA_HCFC1_rep1_mm10_statistics.txt",
#               "db/count_tables/PROseq/HCFC1/stats/AID-Hcfc1-cl17_3hrIAA_HCFC1_rep1_mm10_statistics.txt")
# spikeStats <- c("db/count_tables/PROseq/HCFC1/stats/AID-Hcfc1-cl4_0hrIAA_HCFC1_rep1_dm3_spikeIn_statistics.txt",
#                 "db/count_tables/PROseq/HCFC1/stats/AID-Hcfc1-cl4_3hrIAA_HCFC1_rep1_dm3_spikeIn_statistics.txt",
#                 "db/count_tables/PROseq/HCFC1/stats/AID-Hcfc1-cl17_0hrIAA_HCFC1_rep1_dm3_spikeIn_statistics.txt",
#                 "db/count_tables/PROseq/HCFC1/stats/AID-Hcfc1-cl17_3hrIAA_HCFC1_rep1_dm3_spikeIn_statistics.txt")
# names <- c("control_rep1",
#            "IAA_3h_rep1",
#            "control_rep2",
#            "IAA_3h_rep2")
# conditions <- c("control",
#                 "IAA_3h",
#                 "control",
#                 "IAA_3h")
# controls <- c("control",
#               "control",
#               "control",
#               "control")
# ddsPrefix <- "db/dds/PROseq/HCFC1/transcript"
# fcOutDir <- "db/FC_tables/PROseq/HCFC1/transcript/"

# Parse variables ----
counts <- unlist(tstrsplit(args[1], ","))
refStats <- unlist(tstrsplit(args[2], ","))
spikeStats <- unlist(tstrsplit(args[3], ","))
names <- unlist(tstrsplit(args[4], ","))
names <- gsub("-", ".", names) # Names and conditions do not tolerate "-"
conditions <- unlist(tstrsplit(args[5], ","))
conditions <- gsub("-", ".", conditions) # Names and conditions do not tolerate "-"
controls <- unlist(tstrsplit(args[6], ","))
controls <- gsub("-", ".", controls) # Names and conditions do not tolerate "-"
ddsPrefix <- args[7]
ddsStats <- args[8]
fcOutDir <- args[9]
MAprefix <- args[10]

# Import data ----
dat <- lapply(counts, fread)
names(dat) <- names
dat <- rbindlist(dat, idcol = "condition")
dat[, condition:= factor(condition, unique(condition))]
DF <- dcast(dat, ID~condition, value.var = "count")
DF <- data.frame(DF[, -1], row.names = DF$ID)

# Remove low count reads ----
DF <- DF[rowSums(DF >= 10) >= 2,]

# SampleTable ----
sampleTable <- data.frame(condition = conditions,
                          row.names = names)

# Import read counts and compute normalization factors ----
# Ref genome stats
ref <- lapply(refStats, fread)
names(ref) <- names
ref <- rbindlist(ref, idcol = "sample")
# Spike stats
spike <- lapply(spikeStats, fread)
names(spike) <- names
spike <- rbindlist(spike, idcol = "sample")
# Merge and compute scaling factors
stats <- merge(ref, spike, by= "sample", suffixes= c("", "_spikeIn"))
stats[, cdition:= conditions]
stats[, ctl:= controls]
stats[, spikeInPerc:= round(umi_counts_spikeIn/umi_counts*100, 1)]
stats[, libSizeNorm:= umi_counts/median(umi_counts)]
stats[, spikeInNorm:= umi_counts_spikeIn/median(umi_counts_spikeIn)]

# Vanja's method can only be used with a unique control condition (she was creating one dds object for each pairwise comparison) ----
# Here, I compute it using unique control/condition combination, but I don't know how this would affect objects with multiple combinations
if(all(stats[, length(unique(ctl)), cdition]$V1==1))
{
  stats[, SIctl:= stats[.BY, sum(umi_counts_spikeIn)/sum(umi_counts), on= "cdition==ctl"], ctl]
  stats[, SIsample:= stats[.BY, sum(umi_counts_spikeIn)/sum(umi_counts), on= "cdition"], cdition]
  stats[, combinedNorm:= (SIctl/SIsample-SIctl)/(1-SIctl)]
  stats[, combinedNorm:= combinedNorm*1/umi_counts*10e6]
  stats[, combinedNorm:= max(combinedNorm)/combinedNorm]
}

# DESeq2 analysis ----
# Possible normalization values are "default" (DESeq2 default), "libSize", "spikeIn", "combined" (Vanja's method: note that changing the control sample will change the outcome)
for(normalization in c("libSize", "spikeIn"))
{
  print(paste("Start", normalization, "normalization"))
  # Create DESeq2 object ----
  dds <- DESeqDataSetFromMatrix(countData = DF,
                                colData = sampleTable,
                                design = ~ condition)
  # SizeFactors ----
  if(normalization=="libSize")
    sizeFactors(dds) <- stats$libSizeNorm
  if(normalization=="spikeIn")
    sizeFactors(dds) <- stats$spikeInNorm
  if(normalization=="combined")
    if("combinedNorm" %in% names(stats))
      sizeFactors(dds) <- stats$combinedNorm else
        warning("Provided cdition/control combinations are not unique -> Vanja's combined normalization will be skipped.")
  
  # Compute model and save object ----
  dds <- DESeq(dds)
  saveRDS(dds,
          paste0(ddsPrefix, "_DESeq2_", normalization, "_norm.dds"))
  
  # Import combinations ----
  cmb <- data.table(V1= conditions,
                    V2= controls)
  cmb <- unique(cmb[V1!=V2])
  cmb[, {
    res <- results(dds,
                   contrast = c("condition", V1, V2))
    res <- as.data.frame(res)
    res <- as.data.table(res, keep.rownames = T)
    res[, diff:= fcase(padj<0.05 & log2FoldChange>log2(1.5), "Up-regulated",
                       padj<0.05 & log2FoldChange<(-log2(1.5)), "Down-regulated",
                       default = "Unaffected")]
    outputFile <- paste0(fcOutDir, V1, "_vs_", V2, "_DESeq2_", normalization, "_norm.txt")
    fwrite(res,
           outputFile,
           sep="\t",
           na = NA)
    print(paste(outputFile, "saved"))
    sumup <- table(res$diff)
    print(paste(names(sumup), "= ", sumup, collapse = " | "))
    
    # MA plot ----
    outputPdf <- paste0(MAprefix, "_", V1, "_vs_", V2, "_MA_plot_",  normalization, "_norm.pdf")
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
              main= paste(V1, "vs.", V2, "\n", normalization, " norm."))
    dev.off()
    print(paste(outputPdf, "saved"))
  }, .(V1, V2)]
}

# Plot reads statistics ----
Cc <- c("grey80", "grey60", "grey40", "grey20")
pdf(ddsStats, width = 9)
par(las= 1,
    mai= c(.9,5,.9,2),
    tcl= -0.1,
    mgp= c(1, 0.35, 0))
stats[, {
  bar <- barplot(total,
                 horiz= T,
                 names.arg = paste0(sample, " (", spikeInPerc, "% spike-in)"),
                 xlab= "Number of reads",
                 col= Cc[1])
  barplot(mapped,
          horiz= T,
          add= T,
          col= Cc[2],
          xaxt= "n",
          yaxt= "n")
  barplot(umi_counts,
          horiz= T,
          add= T,
          col= Cc[3],
          xaxt= "n",
          yaxt= "n")
  barplot(umi_counts_spikeIn,
          horiz= T,
          add= T,
          col= Cc[4],
          xaxt= "n",
          yaxt= "n")
  legend(par("usr")[2],
         par("usr")[4],
         legend = c("Total counts", "mapped counts", "UMI counts", "Spike-in"),
         fill= Cc,
         xpd= T,
         bty= "n")
}]
dev.off()
print(paste(ddsStats, "printed! Done!"))
