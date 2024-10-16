#!/usr/bin/env Rscript
suppressMessages(library(vlfunctions, warn.conflicts = FALSE))
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))

# Check arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript script.R {replicates .narrowPeak files (comma separated)} {merged .narrowPeak file} {output_file.narrowPeak}")
}

# args <- c("/groups/stark/vloubiere/projects/epigenetic_cancer/db/peaks/ATAC/ATAC_PH18_rep1_peaks.narrowPeak,/groups/stark/vloubiere/projects/epigenetic_cancer/db/peaks/ATAC/ATAC_PH18_rep2_peaks.narrowPeak",
#           "/groups/stark/vloubiere/projects/epigenetic_cancer/db/peaks/ATAC/ATAC_PH18_merge_peaks.narrowPeak",
#           "db/peaks/CUTNRUN/confident/ATAC_PH18_merge_confident_peaks.narrowPeak")

importNarrowpeak <- function(x)
{
  # Step 1: Import the file
  gr <- rtracklayer::import(x)
  
  # Step 2: Check if the imported object is empty
  if (length(gr) == 0) {
    # Step 3: Create an empty GRanges object with the correct metadata columns
    gr <- GRanges(seqnames = character(),
                  ranges = IRanges(start = integer(), end = integer()),
                  strand = character(),
                  signalValue = numeric(),
                  pValue = numeric(),
                  qValue = numeric(),
                  peak = integer())
  }
  return(as.data.table(gr))
}

# Import files and parse arguments ----
reps <- strsplit(args[[1]], ",")[[1]]
reps <- lapply(reps, function(x) importNarrowpeak(x))
merge <- importNarrowpeak(args[[2]])
conf <- data.table::copy(merge)
output <- args[[3]]

# Overlap replicates ----
merge[, ID:= .I]
ov <- lapply(reps, function(x) vl_intersectBed(merge, x)[, ID])
names(ov) <- paste("Rep", seq(ov))

# Print overlap ----
pdfDir <- "pdf/CUTNRUN/overlap_replicates"
dir.create(pdfDir, showWarnings = F, recursive = T)
pdf(paste0(pdfDir, "/", gsub(".narrowPeak$", ".pdf", basename(output))), 4, 4)
vl_par(mai= c(1.9, 1.9, .9, .9),
       las= 1,
       tcl= -0.1,
       mgp= c(1.5, 0.35, 0),
       cex= 1,
       cex.lab= 9/12,
       cex.axis= 7/12,
       bty= "n",
       lend= 2)
if(sum(lengths(ov))==0)
{
  plot.new()
  text(.5, .5, "No peaks found!")
}else
  vl_upset_plot(ov)
dev.off()

# Extract confident peaks found in all replicates and save ----
for(i in seq(reps))
  conf <- vl_intersectBed(conf, reps[[i]])

# Save ----
if(nrow(conf))
{
  conf[, start:= start-1] # 0-based bed
  conf <- conf[, .(seqnames, start= start-1, end, name, score, strand, signalValue, pValue, qValue, peak)]
}
fwrite(conf,
       output,
       col.names = F,
       quote= F,
       sep= "\t",
       na = ".")
