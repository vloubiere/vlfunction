#' @param obj An object of class vl_enr
#'
#' @param padj.cutoff padjust cutoff to be applied before plotting
#' @param top.enrich Top enrichments to plot (based on padj). Default is to show all enriched features
#' @param order Value to be used for ordering before selecting top enriched. Possible values are "padj", "log2OR". defaut= "padj"
#' @param xlab Default to "Odd Ratio (log2)"
#' @param col Color scale
#' @param breaks Color breaks to be used. Defaults to range of filtered padj.
#'
#' @describeIn vl_motif_enrich method to plot enrichment objects (containing variable, log2OR and padj)
#' @export
plot.vl_enr <- function(obj,
                        padj.cutoff= 0.05,
                        top.enrich= Inf,
                        order= "padj",
                        xlab= "Odd Ratio (log2)",
                        col= c("blue", "red"),
                        breaks= NULL)
{
  if(!(order %in% c("padj", "log2OR")))
    stop("Possible values for order are 'padj', 'log2OR'")
  DT <- data.table::copy(obj)
  # Handle infinite
  if(any(is.infinite(DT$log2OR)))
    warning("Attempt to cap infinite log2OR values max/min finite log2OR. If plot fails, try another representation or cap manually")
  if(any(DT$log2OR==Inf) && nrow(DT[log2OR>0 & is.finite(log2OR)]))
    DT[log2OR==Inf, log2OR:= max(DT[log2OR>0 & is.finite(log2OR), log2OR])]
  if(any(DT$log2OR==(-Inf)) && nrow(DT[log2OR<0 & is.finite(log2OR)]))
    DT[log2OR==(-Inf), log2OR:= min(DT[log2OR<0 & is.finite(log2OR), log2OR])]
  # padj cutoff
  DT <- DT[padj<=padj.cutoff]
  if(nrow(DT))
  {
    # Order
    if(order=="padj")
      setorderv(DT, "padj") else if(order=="log2OR")
        DT <- DT[order(-abs(log2OR))]
    # select top.enrich
    if(nrow(DT)>top.enrich)
      DT <- DT[seq(nrow(DT))<=top.enrich]
    # Plot
    if(is.null(breaks))
    {
      breaks <- range(-log10(DT$padj), na.rm= T)
      if(length(unique(breaks))==1)
        breaks <- breaks+c(-0.5,0.5)
      breaks <- seq(min(breaks), max(breaks), length.out= length(col))
    }
    Cc <- circlize::colorRamp2(breaks, col)
    # If order= 'padj', reorder by Log2OR before plotting
    if(order=="padj")
      setorderv(DT, "log2OR")
    # Barplot
    DT[, y:= barplot(log2OR,
                     horiz= T,
                     names.arg= name,
                     border= NA,
                     col= Cc(-log10(padj)),
                     las= 1,
                     xlab= xlab)]
    # Plot heatkey
    vl_heatkey(breaks = breaks,
               top = DT[.N, y],
               left = par("usr")[2]+strwidth("M"),
               col = col,
               main = "FDR (-log10)")
  }else
    warning("No enrichment found with current cutoffs!")
  # Return
  invisible(DT)
}

#' @describeIn vl_motif_cl_enrich method to plot enrichment objects (containing variable, log2OR and padj)
#' @export
plot.vl_enr_cl <- function(obj,
                           x.breaks,
                           padj.cutoff= 0.05,
                           top.enrich= Inf,
                           order= "padj",
                           color.breaks,
                           cex.balloons= 1,
                           col= c("blue", "red"),
                           main= NA,
                           plot.empty.clusters= T)
{
  if(!(order %in% c("padj", "log2OR")))
    stop("Possible values for order are 'padj', 'log2OR'")
  DT <- data.table::copy(obj)
  # Handle infinite (negative are irrelevant with this plot)
  if(any(is.infinite(DT$log2OR)))
    warning("Attempt to cap infinite log2OR values max/min finite log2OR. If plot fails, try another representation or cap manually")
  if(any(DT$log2OR==Inf) && nrow(DT[log2OR>0 & is.finite(log2OR)]))
    DT[log2OR==Inf, log2OR:= max(DT[log2OR>0 & is.finite(log2OR), log2OR])]
  # Apply cutoffs
  DT <- DT[padj <= padj.cutoff & log2OR > 0]
  if(nrow(DT))
  {
    # Order
    if(order=="padj")
      setorderv(DT, c("cl", "padj")) else if(order=="log2OR")
        DT <- DT[order(cl, -abs(log2OR))]
    # select top.enrich
    if(any(DT[, .N, cl]$N>top.enrich))
      DT <- DT[variable %in% DT[rowid(DT$cl)<=top.enrich, variable]] 
    # Save ordering before dcast
    DT[, variable:= factor(variable, levels= unique(variable))]
    # Remove empty clusters
    if(!plot.empty.clusters)
      DT[, cl:= droplevels(cl)]
    # Add y coordinates to DT and return (matrix upside down)
    DT[, y:= max(as.numeric(variable))-as.numeric(variable)+1]
    # dcast 
    x <- dcast(DT, variable~cl, value.var = "log2OR", drop= F)
    x <- as.matrix(x, 1)
    rownames(x) <- DT[rownames(x), name, on= "variable", mult= "first"]
    color.var <- dcast(DT, variable~cl, value.var = "padj", drop= F)
    color.var <- as.matrix(color.var, 1)
    color.var <- -log10(color.var)
    # Plot
    vl_balloons_plot(x= x,
                     color.var= color.var,
                     x.breaks= x.breaks,
                     col= col,
                     cex.balloons= cex.balloons,
                     main= main,
                     balloon.size.legend= "OR (log2)",
                     balloon.col.legend= "padj (-log10)")
  }else
    warning("No enrichment found with current cutoffs!")
  
  invisible(DT)
}
