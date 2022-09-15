#' @param obj An object of class vl_enr
#'
#' @param padj_cutoff padjust cutoff to be applied before plotting
#' @param top_enrich Top enrichments to plot (based on padj)
#' @param xlab Default to "Odd Ratio (log2)"
#' @param col Color scale
#' @param ... Extra args to be passed to barplot
#'
#' @describeIn vl_motif_enrich method to plot enrichment objects (containing variable, log2OR and padj)
#' @export
plot.vl_enr <- function(obj,
                        padj_cutoff= 0.05,
                        top_enrich= Inf,
                        xlab= "Odd Ratio (log2)",
                        col= c("blue", "red"),
                        ...)
{
  DT <- data.table::copy(obj)
  # Handle infinite
  if(any(!is.finite(DT$log2OR)))
  {
    message("Non finite log2OR values capped to max finite log2OR")
    DT[log2OR==Inf, log2OR:= max(DT[is.finite(log2OR), log2OR])]
    DT[log2OR==(-Inf), log2OR:= min(DT[is.finite(log2OR), log2OR])]
  }
  # padj cutoff
  DT <- na.omit(DT[padj<=padj_cutoff])
  if(nrow(DT)==0)
    stop("No enrichment found with current cutoffs!")
  # select top_enrich
  setorderv(DT, "padj")
  DT <- DT[seq(nrow(DT))<=top_enrich]
  # Plot
  breaks <- range(-log10(DT$padj), na.rm= T)
  Cc <- circlize::colorRamp2(breaks, col)
  setorderv(DT, "log2OR")
  DT[, bar:= barplot(log2OR,
                     horiz= T,
                     names= variable,
                     border= NA,
                     col= Cc(-log10(padj)),
                     las= 1,
                     xlab= xlab,
                     ...)]
  vl_heatkey(breaks = breaks,
             top = last(DT$bar),
             left = par("usr")[2]+strwidth("M"),
             col = col,
             main = "FDR (-log10)")
  invisible(DT)
}

#' @describeIn vl_motif_cl_enrich method to plot cluster enrichment objects (containing variable, log2OR and padj)
#' @export
plot.vl_enr_cl <- function(obj,
                           x_breaks,
                           padj_cutoff= 0.05,
                           top_enrich= Inf,
                           color_breaks,
                           cex.balloons= 1,
                           col= c("cornflowerblue", "lightgrey", "tomato"),
                           main= NA)
{
  DT <- data.table::copy(obj)
  # Handle infinite
  if(any(!is.finite(DT$log2OR)))
  {
    message("Non finite log2OR values capped to max finite log2OR")
    DT[log2OR==Inf, log2OR:= max(DT[is.finite(log2OR), log2OR])]
    DT[log2OR==(-Inf), log2OR:= min(DT[is.finite(log2OR), log2OR])]
  }
  # Apply cutoffs
  DT <- DT[padj <= padj_cutoff & log2OR > 0]
  if(nrow(DT)==0)
    stop("No enrichment found with current cutoffs!")
  # select top_enrich
  setorderv(DT, c("cl", "padj"))
  sel <- DT[rowid(DT$cl)<=top_enrich, variable]
  DT <- DT[variable %in% sel]
  # dcast before plotting
  DT[, variable:= factor(variable, levels= unique(variable))]
  DT[, cl:= droplevels(cl)]
  x <- dcast(DT, variable~cl, value.var = "log2OR", drop= F)
  x <- as.matrix(x, 1)
  color_var <- dcast(DT, variable~cl, value.var = "padj", drop= F)
  color_var <- as.matrix(color_var, 1)
  color_var <- -log10(color_var)
  # Plot
  pl <- match.call()
  pl$obj <- pl$padj_cutoff <- pl$top_enrich <- NULL
  pl$x <- x
  pl$color_var <- color_var
  pl$cex.balloons <- cex.balloons
  pl$balloon_size_legend <- "OR (log2)"
  pl$balloon_col_legend <- "padj (-log10)"
  pl[[1]] <- quote(vl_balloons_plot)
  eval(pl, envir = parent.frame())
  
  invisible(list(x= x,
                 color_var= color_var))
}