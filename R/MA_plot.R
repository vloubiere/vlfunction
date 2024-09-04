#' Plots a nice MA plot
#'
#' @param DESeq2_dat A data.table containing a Deseq2 output sheet with log2FoldChange and baseMean.
#' @param padj.cutoff padjust cutoff to be used to call up/down genes. Default= 0.05.
#' @param log2FC.cutoff log2FoldChange cutoff to be used to call up/down genes. Default= log2(1.5).
#' @param main Title
#'
#' @return MA plot
#' @export
#'
#' @examples
#' vl_MAplot(fread("/groups/stark/vloubiere/projects/epigenetic_cancer/db/FC_tables/RNA/epiCancer_GFP_PH29_vs_W29.txt"),
#'           main= "PH29")
vl_MAplot <- function(DESeq2_dat, main= NA, padj.cutoff= 0.05, log2FC.cutoff= log2(1.5))
{
  # Compute diff column
  if("diff" %in% DESeq2_dat)
    DESeq2_dat$diff <- NULL
  DESeq2_dat[log2FoldChange>log2FC.cutoff & padj<padj.cutoff, diff:= "Up-regulated"]
  DESeq2_dat[log2FoldChange<(-log2FC.cutoff) & padj<padj.cutoff, diff:= "Down-regulated"]
  
  # Compute diff column
  DESeq2_dat[, {
    clip <- quantile(log2FoldChange, c(0.001, 0.999))
    y <- log2FoldChange
    y[y<clip[1]] <- clip[1]
    y[y>clip[2]] <- clip[2]
    Cc <- fcase(diff=="Up-regulated", "tomato",
                diff=="Down-regulated", "cornflowerblue",
                default = "lightgrey")
    nUp <- sum(diff=="Up-regulated")
    nUp <- formatC(nUp, big.mark = ",")
    nDown <- sum(diff=="Down-regulated")
    nDown <- formatC(nDown, big.mark = ",")
    plot(log10(baseMean),
         y,
         col= adjustcolor(Cc, .5),
         pch= ifelse(y!=log2FoldChange, 17, 16),
         ylab= "PRO-Seq fold change (log2)",
         frame= F,
         xaxt= "n",
         cex= .5,
         main= main)
    axis(1, padj= -1.45)
    abline(h= 0, lty= 3)
    legend(par("usr")[2],
           par("usr")[4],
           col= adjustcolor(c("tomato", "cornflowerblue"), 0.5),
           legend= c(paste0("Up-regulated (", nUp, ")"), 
                     paste0("Down-regulated (", nDown, ")")),
           pch= 16,
           bty= "n",
           cex= 6/12,
           border= NA,
           xpd= NA)
  }]
}
