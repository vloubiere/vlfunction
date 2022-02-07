#' Counts motifs
#'
#' Counts motif occurences in a set of regions
#'
#' @param sequences Sequences to analyse. If provided takes over bed argument (in the case where both are specified)
#' @param bed Either a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param genome Genome to be used for coordinates ("dm6, "dm3")
#' @param sel List of motifs to compute. see vl_Dmel_motifs_DB_full$motif
#' 
#' @examples 
#' test <- cl_motif_counts(bed = GRanges("chr3R", IRanges(c(2e6, 3e6), c(2e6, 3e6)+1e3)),
#' genome= "dm6",
#' resize= T, 
#' extend = c(500, 500),
#' sel= vl_Dmel_motifs_DB$metadata$motif_name[1:10])
#' 
#' @return Network plot.
#' @export

vl_motif_counts <- function(sequences= NULL,
                            bed= NULL, 
                            genome,
                            sel= vl_Dmel_motifs_DB_full[!is.na(vl_Dmel_motifs_DB_full$FBgn), motif])
{
  # Checks
  if(is.null(sequences))
  {
    if(is.null(bed))
      stop("sequences or bed+genome must be specified") else
    {
      if(!vl_isDTranges(bed))
        bed <- vl_importBed(bed)
      sequences <- BSgenome::getSeq(BSgenome::getBSgenome(genome), bed$seqnames, bed$start, bed$end, as.character= T)
    }
  } 
  
  # Select motifs
  if(any(!sel %in% vl_Dmel_motifs_DB_full$motif))
    stop("Some motif provided in 'sel' do not exist in vl_Dmel_motifs_DB_full$motif")
  sub <- vl_Dmel_motifs_DB_full[motif %in% sel]
  mot <- do.call(TFBSTools::PWMatrixList, 
                 sub$pwms_log_odds)
  
  # Count motifs
  res <- as.data.table(as.matrix(motifmatchr::matchMotifs(mot,
                                                          sequences,
                                                          p.cutoff= 5e-4,
                                                          bg= "even",
                                                          out= "scores")@assays@data[["motifCounts"]]))
  res <- as.data.table(res)
  setnames(res, sub$motif)
  return(res)  
}


#' Compute motif enrichment
#'
#' Compute motif enrichment for the cluster in cl_columns, using all the lines as background
#'
#' @param obj A fisher test matrix similar to ?vl_motif_cl_enrich() ouput.
#' @param x_breaks Breaks to be used for balloons'sizes
#' @param padj_cutoff cutoff for ballons to be ploted
#' @param log2OR_cutoff cutoff for ballons to be ploted
#' @param N_top Select top enriched motifs/cluster
#' @param color_breaks Color breaks used for coloring
#' @param col Vector of colors used for coloring
#' @param main Title. Default= NA
#' @param auto_margins Use auto margins? Default= T
#' @return ballons plot with padj and log2OR
#' @export

vl_motif_cl_enrich_plot_only <- function(obj,
                                         x_breaks,
                                         padj_cutoff= 0.05,
                                         log2OR_cutoff= 0,
                                         N_top= Inf,
                                         color_breaks,
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
  DT <- DT[, .(cl, log2OR, padj, any(padj <= padj_cutoff & log2OR > log2OR_cutoff)), motif][(V4), !"V4"]
  if(nrow(DT)==0)
    stop("No enrichment found with provided padj cutoff!")
  # Format plotting object
  DT[, '-log10(padj)':= -log10(padj)]
  # Order and select top motif/cluster
  setorderv(DT, 
            c("cl", "-log10(padj)", "log2OR", "motif"), 
            order = c(1, -1, -1, 1))
  DT[log2OR>0, rank:= seq(.N), cl]
  DT <- DT[motif %in% DT[rank<=N_top, motif]]
  # Handle infinite log2OR
  if(any(is.infinite(DT$log2OR)))
  {
    warning("There are suspcious infinite log2OR values found after padj cutoff -> capped to max(finite)")
    print(unique(DT[is.infinite(DT$log2OR), motif]))
    DT[log2OR==Inf, log2OR:= max(DT[is.finite(log2OR), log2OR])]
    DT[log2OR==(-Inf), log2OR:= min(DT[is.finite(log2OR), log2OR])]
  }
  # Remove values that were useful for earlier ordering but do not meet minimum criteria
  DT <- DT[padj<0.05 & log2OR>=0]
  DT[, motif:= factor(motif, levels= unique(DT$motif))]
  
  #-----------------------------#
  # PLOT
  #-----------------------------#
  x <- as.matrix(dcast(DT, motif~cl, value.var = "log2OR"), 1)
  color_var <- as.matrix(dcast(DT, motif~cl, value.var = "-log10(padj)"), 1)
  
  pl <- match.call()
  pl$obj <- pl$padj_cutoff <- pl$log2OR_cutoff <- pl$N_top <- NULL
  pl$x <- x
  pl$color_var <- color_var
  pl$balloon_size_legend <- "OR (log2)"
  pl$balloon_col_legend <- "padj (-log10)"
  pl[[1]] <- quote(vl_balloons_plot)
  eval(pl)
}

#' Compute motif enrichment
#'
#' Compute motif enrichment for the cluster in cl_columns, using all the lines as background
#'
#' @param counts_matrix A matrix containing motif counts (rows= sequences, cols= motifs)
#' @param cl_IDs vector of cluster IDs (used to split the data.table)
#' @param control_cl IDs of clusters to be used as background. default to unique(cl_IDs), meaning all clusters are used except the one being tested
#' @param plot Should the result be plot using balloons plot?
#' @param padj_cutoff cutoff for ballons to be ploted
#' @param log2OR_cutoff cutoff for ballons to be ploted
#' @param N_top Select top enriched motifs/cluster
#' @param x_breaks Breaks used for ballon's sizes
#' @param color_breaks Color breaks used for coloring
#' @param col Vector of colors used for coloring
#' @param main Title. Default= NA
#' @param auto_margins Use auto margins? Default= T
#' @return Fisher test data.table.
#' @export

vl_motif_cl_enrich <- function(counts_matrix, 
                               cl_IDs,
                               control_cl= unique(cl_IDs),
                               plot= F,
                               padj_cutoff= 0.00001,
                               log2OR_cutoff= 0,
                               N_top= Inf,
                               x_breaks,
                               color_breaks,
                               col= c("cornflowerblue", "lightgrey", "tomato"),
                               main= NA,
                               auto_margins = T)
{
  if(!is.matrix(counts_matrix))
    stop("counts_matrix should be a matrix")
  if(any(sapply(seq(ncol(counts_matrix)), function(i) !is.numeric(counts_matrix[,i]))))
    stop("counts_matrix should only contain numeric values")
  if(is.null(rownames(counts_matrix)))
    rownames(counts_matrix) <- seq(nrow(counts_matrix))
  counts_matrix <- as.data.table(counts_matrix, keep.rownames= T)
  
  # Format table
  res <- melt(cbind(cl= cl_IDs, counts_matrix), 
              id.vars = c("rn", "cl"),
              variable.name = "motif")
  # Enrichment
  res <- res[, {
    counts <- data.table(ccl= cl, 
                         cval= value,
                         key= "ccl")
    .c <- .SD[, {
      # Contingency table restricted to control cluster(s) and the tested one
      tab <- table(counts[.(unique(c(cl, control_cl))), .(ccl==cl, cval>0)])
      if(identical(c(2L,2L), dim(tab))) # Filter out motif for which all counts fall within one category
        fisher.test(tab)[c("estimate", "p.value")] else
          list(estimate= as.numeric(NA), `p.value`= as.numeric(NA))
    }, cl]
  }, motif]
  setnames(res, c("estimate", "p.value"), c("OR", "pval"))
  # padj...
  res[, padj:= p.adjust(pval, method = "fdr"), pval]
  res[, log2OR:= log2(OR)]
  res <- res[, .(motif, cl, log2OR, padj)]
  # Plot
  if(plot)
  {
    pl <- match.call()
    pl$counts_matrix <- pl$cl_IDs <- pl$control_cl <- pl$plot <- NULL
    pl$obj <- res
    pl[[1]] <- quote(vl_motif_cl_enrich_plot_only)
    eval(pl)
  }
  
  return(res)
}
