MA_plot <- function(DESeq2_dat, pdf, main)
{
  pdf(pdf, 4, 3)
  par(mai= c(.9,1.5,.9,1.3),
      cex.axis= 6/12,
      cex.lab= 7/12,
      las= 1,
      tcl= -0.1,
      mgp= c(.8, 0.25, 0),
      font.main= 1,
      cex.main= 9/12)
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
         cex= .5)
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
  title(main= main)
  dev.off()
}