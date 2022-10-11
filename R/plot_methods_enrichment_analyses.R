#' @param obj An object of class vl_enr
#'
#' @param padj_cutoff padjust cutoff to be applied before plotting
#' @param top_enrich Top enrichments to plot (based on padj)
#' @param order Value to be used for ordering before selecting top enriched. Possible values are "padj", "log2OR". defaut= "padj"
#' @param xlab Default to "Odd Ratio (log2)"
#' @param col Color scale
#' @param breaks Color breaks to be used. Defaults to range of filtered padj.
#'
#' @describeIn vl_motif_enrich method to plot enrichment objects (containing variable, log2OR and padj)
#' @export
plot.vl_enr <- function(obj,
                        padj_cutoff= 0.05,
                        top_enrich= NA,
                        order= "padj",
                        xlab= "Odd Ratio (log2)",
                        col= c("blue", "red"),
                        breaks= NULL)
{
  if(!(order %in% c("padj", "log2OR")))
    stop("Possible values for order are 'padj', 'log2OR'")
  DT <- data.table::copy(obj)
  # Handle infinite
  if(any(!is.finite(DT$log2OR)))
    warning("Non finite log2OR values capped to max finite log2OR")
  if(any(DT$log2OR==Inf))
    if(any(is.finite(DT[log2OR>0, log2OR])))
      DT[log2OR==Inf, log2OR:= max(DT[log2OR>0 & is.finite(log2OR), log2OR])] else
        stop("Inf OR and no finite pos OR to use for capping. Use other visualization!")
  if(any(DT$log2OR==(-Inf)))
    if(any(is.finite(DT[log2OR<0, log2OR])))
      DT[log2OR==(-Inf), log2OR:= min(DT[log2OR<0 & is.finite(log2OR), log2OR])] else
        stop("-Inf OR and no finite neg OR to use for capping. Use other visualization!")
  # padj cutoff
  DT <- DT[padj<=padj_cutoff]
  if(nrow(DT)>0)
  {
    # select top_enrich
    if(!is.na(top_enrich) && nrow(DT)>top_enrich)
    {
      if(order=="padj")
        setorderv(DT, "padj") else if(order=="log2OR")
          DT <- DT[order(-abs(log2OR))]
      DT <- DT[seq(nrow(DT))<=top_enrich]
    }
    # Plot
    if(is.null(breaks))
    {
      breaks <- range(-log10(DT$padj), na.rm= T)
      if(length(unique(breaks))==1)
        breaks <- breaks+c(-0.5,0.5)
      breaks <- seq(min(breaks), max(breaks), length.out= length(col))
    }
    Cc <- circlize::colorRamp2(breaks, col)
    # Barplot
    setorderv(DT, "log2OR")
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
                           x_breaks,
                           padj_cutoff= 0.05,
                           top_enrich= NA,
                           order= "padj",
                           color_breaks,
                           cex.balloons= 1,
                           col= c("blue", "red"),
                           main= NA,
                           plot_empty_clusters= T)
{
  if(!(order %in% c("padj", "log2OR")))
    stop("Possible values for order are 'padj', 'log2OR'")
  DT <- data.table::copy(obj)
  # Handle infinite
  if(any(!is.finite(DT$log2OR)))
    warning("Non finite log2OR values capped to max finite log2OR")
  if(any(DT$log2OR==Inf))
    if(any(is.finite(DT[log2OR>0, log2OR])))
      DT[log2OR==Inf, log2OR:= max(DT[log2OR>0 & is.finite(log2OR), log2OR])] else
        stop("Inf OR and no finite pos OR to use for capping. Use other visualization!")
  if(any(DT$log2OR==(-Inf)))
    if(any(is.finite(DT[log2OR<0, log2OR])))
      DT[log2OR==(-Inf), log2OR:= min(DT[log2OR<0 & is.finite(log2OR), log2OR])] else
        stop("-Inf OR and no finite neg OR to use for capping. Use other visualization!")
  # Apply cutoffs
  DT <- DT[padj <= padj_cutoff & log2OR > 0]
  if(nrow(DT)>0)
  {
    # select top_enrich
    if(!is.na(top_enrich) && nrow(DT)>top_enrich)
    {
      if(order=="padj")
        setorderv(DT, c("cl", "padj")) else if(order=="log2OR")
          DT <- DT[order(cl, -abs(log2OR))]
        sel <- DT[rowid(DT$cl)<=top_enrich, variable]
        DT <- DT[variable %in% sel]
    }
    # Save ordering before dcast
    setorderv(DT, c("cl", "padj"))
    DT[, variable:= factor(variable, levels= unique(variable))]
    if(!plot_empty_clusters)
      DT[, cl:= droplevels(cl)]
    # Add y coordinates to DT and return
    DT[, y:= as.numeric(variable)]
    # dcast 
    x <- dcast(DT, variable~cl, value.var = "log2OR", drop= F)
    x <- as.matrix(x, 1)
    rownames(x) <- DT[rownames(x), name, on= "variable", mult= "first"]
    color_var <- dcast(DT, variable~cl, value.var = "padj", drop= F)
    color_var <- as.matrix(color_var, 1)
    color_var <- -log10(color_var)
    rownames(color_var) <- DT[rownames(color_var), name, on= "variable", mult= "first"]
    # Plot
    vl_balloons_plot(x= x,
                     color_var= color_var,
                     x_breaks= x_breaks,
                     col= col,
                     cex.balloons= cex.balloons,
                     main= main,
                     balloon_size_legend= "OR (log2)",
                     balloon_col_legend= "padj (-log10)")
  }else
    warning("No enrichment found with current cutoffs!")
  
  invisible(DT)
}
