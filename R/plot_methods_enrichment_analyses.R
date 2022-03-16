#' @describeIn vl_motif_enrich method to plot enrichment objects (containing variable, log2OR and padj)
#' @export
plot.vl_enr <- function(obj,
                        padj_cutoff= 0.05,
                        top_enrich= Inf)
{
  DT <- data.table::copy(obj)
  DT <- na.omit(DT[padj<=padj_cutoff][order(log2OR)])
  DT[, `-log10(padj)`:= -log10(padj)]
  DT[, rank:= order(log2OR, decreasing = T)]
  # Handle infinite
  if(any(is.infinite(DT$log2OR)))
  {
    warning("There are suspcious infinite log2OR values found after padj cutoff -> capped to max(finite)")
    print(unique(DT[is.infinite(DT$log2OR), variable]))
    DT[log2OR==Inf, log2OR:= max(DT[is.finite(log2OR) & log2OR>=0, log2OR])]
    DT[log2OR==(-Inf), log2OR:= min(DT[is.finite(log2OR)  & log2OR<=0, log2OR])]
  }
  DT <- DT[rank<=top_enrich]
  breaks <- range(DT$`-log10(padj)`, na.rm= T)
  Cc <- circlize::colorRamp2(breaks,
                             c("blue", "red"))
  DT[, bar:= barplot(log2OR,
                     horiz= T,
                     names= variable,
                     border= NA,
                     col= Cc(`-log10(padj)`),
                     las= 1,
                     xlab= "Odd Ratio (log2)")]
  vl_heatkey(breaks = breaks,
             top = last(DT$bar),
             left = par("usr")[2]+strwidth("M"),
             col = c("blue", "red"),
             main = "FDR (-log10)")
  invisible(DT)
}

#' @describeIn vl_motif_cl_enrich method to plot cluster enrichment objects (containing variable, log2OR and padj)
#' @export
plot.vl_enr_cl <- function(obj,
                           x_breaks,
                           padj_cutoff= 0.05,
                           log2OR_cutoff= 0,
                           N_top= Inf,
                           color_breaks,
                           cex.balloons= 1,
                           col= c("cornflowerblue", "lightgrey", "tomato"),
                           main= NA,
                           auto_margins = T)
{
  DT <- data.table::copy(obj)
  if(log2OR_cutoff<0)
    stop("log2OR_cutoff should be >= 0")
  
  #----------------------------------#
  # Generate plot table
  #----------------------------------#
  # Cutoffs
  DT[, check:= any(padj <= padj_cutoff & log2OR > log2OR_cutoff), variable]
  DT <- DT[(check), .(cl, log2OR, padj), variable]
  if(nrow(DT)==0)
    stop("No enrichment found with provided padj cutoff!")
  # Format plotting object
  DT[, '-log10(padj)':= -log10(padj)]
  # Order and select top variable/cluster
  setorderv(DT, 
            c("cl", "-log10(padj)", "log2OR", "variable"), 
            order = c(1, -1, -1, 1))
  DT[log2OR>0, rank:= seq(.N), cl]
  DT <- DT[variable %in% DT[rank<=N_top, variable]]
  # Handle infinite log2OR/padj
  if(any(is.infinite(DT[, c(log2OR, `-log10(padj)`)])))
  {
    warning("There are suspcious infinite log2OR/padj values found after padj cutoff -> capped to max(finite)")
    cols <- c("log2OR", "-log10(padj)")
    DT[, (cols):= lapply(.SD, function(x) fcase(!is.finite(x) & x<0, min(x[is.finite(x)]),
                                                !is.finite(x) & x>0, max(x[is.finite(x)]),
                                                is.finite(x), x)), .SDcols= cols]
  }
  # Remove values that were useful for earlier ordering but do not meet minimum criteria
  DT <- DT[padj<0.05 & log2OR>=0]
  DT[, variable:= factor(variable, levels= unique(DT$variable))]
  
  #-----------------------------#
  # PLOT
  #-----------------------------#
  x <- as.matrix(dcast(DT, variable~cl, value.var = "log2OR"), 1)
  color_var <- as.matrix(dcast(DT, variable~cl, value.var = "-log10(padj)"), 1)
  
  pl <- match.call()
  pl$obj <- pl$padj_cutoff <- pl$log2OR_cutoff <- pl$N_top <- NULL
  pl$x <- x
  pl$color_var <- color_var
  pl$cex.balloons <- cex.balloons
  pl$balloon_size_legend <- "OR (log2)"
  pl$balloon_col_legend <- "padj (-log10)"
  pl[[1]] <- quote(vl_balloons_plot)
  eval(pl)
}