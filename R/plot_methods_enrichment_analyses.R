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
  DT[, `-log10(padj)`:= -log10(padj)]
  # padj and min set counts cutoff
  DT <- na.omit(DT[padj<=padj_cutoff])
  # Handle infinite
  if(any(DT$log2OR==Inf))
    if(any(DT$log2OR>0 & DT$log2OR<Inf))
    {
      warning("Inf log2OR values capped to max finite log2OR")
      DT[log2OR==Inf, log2OR:= max(DT[is.finite(log2OR), log2OR])]
    }else
      stop("all(log2OR[log2OR>0]==Inf) -> could not be capped")
  if(any(DT$log2OR==(-Inf)))
    if(any(DT$log2OR<0 & DT$log2OR>(-Inf)))
    {
      warning("-Inf log2OR values capped to min finite log2OR")
      DT[log2OR==(-Inf), log2OR:= min(DT[is.finite(log2OR), log2OR])]
    }else
      stop("all(log2OR[log2OR<0]==(-Inf)) -> could not be capped")
  # select top_enrich
  DT[, rank:= order(`-log10(padj)`, decreasing = T)]
  DT <- DT[rank<=top_enrich]
  # Ploting
  breaks <- range(DT$`-log10(padj)`, na.rm= T)
  Cc <- circlize::colorRamp2(breaks, col)
  DT[, bar:= barplot(log2OR,
                     horiz= T,
                     names= variable,
                     border= NA,
                     col= Cc(`-log10(padj)`),
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
                           log2OR_cutoff= 0,
                           top_enrich= Inf,
                           color_breaks,
                           cex.balloons= 1,
                           col= c("cornflowerblue", "lightgrey", "tomato"),
                           main= NA,
                           auto_margins = T)
{
  DT <- data.table::copy(obj)
  if(log2OR_cutoff<0)
    stop("log2OR_cutoff should be >= 0")
  DT[, `-log10(padj)`:= -log10(padj)]
  # Apply cutoffs
  sel <- DT[padj <= padj_cutoff & log2OR > log2OR_cutoff, variable]
  DT <- na.omit(DT[variable %in% sel & log2OR>0])
  if(nrow(DT)==0)
    stop("No enrichment found with provided padj cutoff!")
  
  #----------------------------------#
  # Generate plot table
  #----------------------------------#
  # Handle infinite
  if(any(DT$log2OR==Inf))
    if(any(DT$log2OR>0 & DT$log2OR<Inf))
    {
      warning("Inf log2OR values capped to max finite log2OR")
      DT[log2OR==Inf, log2OR:= max(DT[is.finite(log2OR), log2OR])]
    }else
      stop("all(log2OR[log2OR>0]==Inf) -> could not be capped")
  # Order and select top variable/cluster
  setorderv(DT, 
            c("cl", "-log10(padj)", "log2OR", "variable"), 
            order = c(1, -1, -1, 1))
  DT[, rank:= rowid(cl)]
  DT <- DT[variable %in% DT[rank<=top_enrich, variable]]
  # Remove values that were useful for earlier ordering but do not meet minimum criteria
  DT <- DT[padj<0.05]
  DT[, variable:= factor(variable, levels= unique(DT$variable))]
  
  #-----------------------------#
  # PLOT
  #-----------------------------#
  x <- dcast(DT, variable~cl, value.var = "log2OR", drop= F)
  x <- as.matrix(x, 1)
  color_var <- dcast(DT, variable~cl, value.var = "-log10(padj)", drop= F)
  color_var <- as.matrix(color_var, 1)
  
  pl <- match.call()
  pl$obj <- pl$padj_cutoff <- pl$log2OR_cutoff <- pl$top_enrich <- NULL
  pl$x <- x
  pl$color_var <- color_var
  pl$cex.balloons <- cex.balloons
  pl$balloon_size_legend <- "OR (log2)"
  pl$balloon_col_legend <- "padj (-log10)"
  pl[[1]] <- quote(vl_balloons_plot)
  eval(pl)
  
  invisible(list(x= x,
                 color_var= color_var))
}