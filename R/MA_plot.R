#' Plots a nice MA plot
#'
#' @param data Either a valid path to a '.txt' file or a data.table object corresponding to a DESeq2 output sheet  containing both 'log2FoldChange' and 'baseMean' columns.
#' @param main Title
#' @param padj.cutoff padjust cutoff to be used to call up/down genes. Default= 0.05.
#' @param log2FC.cutoff log2FoldChange cutoff to be used to call up/down genes. Default= log2(1.5).
#'
#' @return MA plot
#' @export
#'
#' @examples
#' vl_MAplot(data= "/groups/stark/vloubiere/projects/epigenetic_cancer/db/FC_tables/RNA/epiCancer_GFP_PH29_vs_W29.txt",
#'           main= "PH29")
vl_MAplot <- function(data,
                      main= NA,
                      padj.cutoff= 0.05,
                      log2FC.cutoff= log2(1.5))
{
  # Checks
  if(!is.data.table(data)){
    if(!is.character(data) || !grepl(".txt$", data))
      stop("data should either be a  path to a '.txt' file or a data.table") else
        data <- fread(data)
  }

  # Compute diff column ----
  if("diff" %in% names(data))
    data$diff <- NULL
  data[, diff:= {
    fcase(log2FoldChange>log2FC.cutoff & padj<padj.cutoff, "Up-regulated",
          log2FoldChange<(-log2FC.cutoff) & padj<padj.cutoff, "Down-regulated",
          default = "Unaffected")
  }]

  # Colors ----
  data[, col:= {
    fcase(diff=="Up-regulated", "tomato",
          diff=="Down-regulated", "cornflowerblue",
          default = "lightgrey")
  }]

  # Ordering: grey points plotted first ----
  data <- data[order(col=="lightgrey", decreasing = TRUE)]

  # Compute diff column ----
  data[, {
    # Clip
    clip <- quantile(log2FoldChange, c(0.001, 0.999))
    y <- log2FoldChange
    y[y<clip[1]] <- clip[1]
    y[y>clip[2]] <- clip[2]

    # Plot
    plot(log10(baseMean),
         y,
         col= adjustcolor(col, .5),
         pch= ifelse(y!=log2FoldChange, 17, 16),
         ylab= "PRO-Seq fold change (log2)",
         frame= F,
         xaxt= "n",
         cex= .5,
         main= main)
    axis(1, padj= -1.45)
    abline(h= 0, lty= 3)

    # Legend
    nUp <- sum(diff=="Up-regulated")
    nUp <- formatC(nUp, big.mark = ",")
    nDown <- sum(diff=="Down-regulated")
    nDown <- formatC(nDown, big.mark = ",")
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
