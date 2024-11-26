#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("Please specify:\n
       [required] 1/ Gene summary output file from mageck\n
       [required] 2/ logFC cutoff for hits\n
       [required] 3/ FDR cutoff for hits\n
       [required] 4/ pdf output file\n")
}

require(ggplot2)
require(ggrepel)
require(data.table)

# Import data ----
gene_summary <- args[1]
dat <- fread(gene_summary)
logFCcutoff <- args[2]
FDRcutoff <- args[3]
pdf <- args[4]

# Identify hits and save ----
setnames(dat, "id", "ORF")
dat[, hit:= `pos|lfc`>logFCcutoff & `pos|fdr`<FDRcutoff]
dat[, hitFDR5:= `pos|fdr`<1e-5]
fwrite(dat,
       gsub(".gene_summary.txt$", "_FC_MAGeCK.txt", gene_summary),
       na = NA,
       sep= "\t")

# Prepare for plotting ----
pl <- data.table::copy(dat)
pl[, logFDR:= -log10(`pos|fdr`)]
yMax <- quantile(pl$`pos|lfc`, .999, na.rm= T)
pl[, shape:= ifelse(logFDR > yMax, "triangle", "circle")]
pl[logFDR>yMax, logFDR:= yMax]
pl[, col:= fcase(hit, "Hit", hitFDR5, "FDR<1e-5", default = "None")]
plMageck <- ggplot(pl, aes(x = `pos|lfc`, y = logFDR)) +
  geom_point(aes(color = col, shape = shape)) +
  geom_text_repel(data = pl[(col!="None")], 
                  aes(label = ORF, col= col), 
                  max.overlaps = Inf,
                  size = 2) +
  theme_minimal() +
  labs(title = "Mageck method",
       x = "Fold Change (log)",
       y = "FDR (-log10)") +
  scale_color_manual(values = c("Hit" = "red", "FDR<1e-5" = "blue", "None" = "lightgrey"), name = "Significant")+
  ylim(0, yMax) +
  guides(shape = "none")  # Remove shape legend

# Print pdf
pdf(pdf)
plot(plMageck)
dev.off()