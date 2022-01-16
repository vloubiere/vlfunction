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
                            genome= "dm6",
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
      sequences <- BSgenome::getSeq(BSgenome::getBSgenome(genome), seqnames, start, end, as.character= T)
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
#' @param padj_cutoff cutoff for ballons to be ploted
#' @param log2OR_cutoff cutoff for ballons to be ploted
#' @param N_top Select top enriched motifs/cluster
#' @param color_breaks Color breaks used for coloring
#' @param col Vector of colors used for coloring
#' @param main Title. Default= NA
#' @param cex.balloons scaling factor for balloons size
#' @param auto_margins Use auto margins? Default= T
#' @return ballons plot with padj and log2OR
#' @export

vl_motif_cl_enrich_plot_only <- function(obj,
                                         padj_cutoff= 0.00001,
                                         log2OR_cutoff= 0,
                                         N_top= Inf,
                                         color_breaks,
                                         col= c("cornflowerblue", "lightgrey", "tomato"),
                                         main= NA,
                                         cex.balloons = 4, 
                                         auto_margins = T)
{
  DT <- data.table::copy(obj)
  if(log2OR_cutoff<0)
    stop("log2OR_cutoff should be > 0")
  
  #----------------------------------#
  # Generate plot table
  #----------------------------------#
  # Cutoffs
  sel <- DT[, any(padj <= padj_cutoff & log2OR > log2OR_cutoff), motif][, V1]
  if(!any(sel))
    stop("No enrichment found with provided paj cutoff!")
  DT <- DT[(sel & padj < 0.05 & log2OR > 0)]
  # Format plotting object
  DT[, '-log10(pval)':= -log10(padj)]
  # Select top motif/cluster
  setorderv(DT, "-log10(pval)", order = -1)
  top_motifs <- DT[, motif[seq(.N)<=N_top], cl]$V1
  DT <- DT[motif %in% top_motifs]
  # Handle infinite log2OR
  if(any(is.infinite(DT$log2OR)))
  {
    warning("There are suspcious infinite log2OR values found after padj cutoff")
    print(unique(DT[is.infinite(DT$log2OR), motif]))
  }
  DT[log2OR==Inf, log2OR:= max(pl[is.finite(log2OR), log2OR])]
  DT[log2OR==(-Inf), log2OR:= min(pl[is.finite(log2OR), log2OR])]
  # Y ordering
  setorderv(DT, 
            c("cl", "-log10(pval)", "log2OR", "motif"), 
            order = c(1, -1, -1, 1))
  DT[, motif:= factor(motif, levels= unique(DT$motif))]
  
  #-----------------------------#
  # PLOT
  #-----------------------------#
  x <- as.matrix(dcast(DT, motif~cl, value.var = "log2OR"), 1)
  color_var <- as.matrix(dcast(DT, motif~cl, value.var = "-log10(pval)"), 1)
  if(missing(color_breaks))
    color_breaks <- seq(min(color_var, na.rm= T), 
                        max(color_var, na.rm= T), 
                        length.out= length(col))
  vl_balloons_plot(x = x,
                   color_var = color_var,
                   color_breaks = color_breaks,
                   col= col,
                   main= main,
                   cex.balloons = cex.balloons,
                   auto_margins = auto_margins,
                   balloon_size_legend= "OR (log2)",
                   balloon_col_legend = "padj (-log10)")
}

#' Compute motif enrichment
#'
#' Compute motif enrichment for the cluster in cl_columns, using all the lines as background
#'
#' @param counts_matrix A matrix containing motif counts (rows= sequences, cols= motifs)
#' @param cl_IDs vector of cluster IDs (used to split the data.table)
#' @param plot Should the result be plot using balloons plot?
#' @param padj_cutoff cutoff for ballons to be ploted
#' @param log2OR_cutoff cutoff for ballons to be ploted
#' @param N_top Select top enriched motifs/cluster
#' @param color_breaks Color breaks used for coloring
#' @param col Vector of colors used for coloring
#' @param main Title. Default= NA
#' @param cex.balloons scaling factor for balloons size
#' @param auto_margins Use auto margins? Default= T
#' @return Fisher test data.table.
#' @export

vl_motif_cl_enrich <- function(counts_matrix, 
                               cl_IDs,
                               plot= F,
                               padj_cutoff= 0.00001,
                               log2OR_cutoff= 0,
                               N_top= Inf,
                               color_breaks,
                               col= c("cornflowerblue", "lightgrey", "tomato"),
                               main= NA,
                               cex.balloons = 4, 
                               auto_margins = T)
{
  if(!is.data.table(counts_matrix))
    counts_matrix <- as.data.table(counts_matrix, keep.rownames= T)
  counts_matrix <- data.table::copy(counts_matrix)
  
  # Format table
  res <- melt(cbind(cl= cl_IDs, counts_matrix), 
              measure.vars = names(counts_matrix), 
              variable.name = "motif")
  res[, value:= value>0]
  res <- merge(dcast(res, motif~value, fun.aggregate = length), 
               dcast(res, motif+cl~value, fun.aggregate = length), 
               suffixes= c("_all", "_cluster"))
  # Compute only sequences for which counts>0
  check <- res[, FALSE_cluster>0 & TRUE_cluster>0]
  res[(check), c("OR", "pval"):= {
    mat <- matrix(unlist(.BY), nrow= 2, byrow = T)
    fisher.test(mat)[c("estimate", "p.value")]
  }, .(TRUE_cluster, FALSE_cluster, TRUE_all, FALSE_all)]
  res[, padj:= p.adjust(pval, method = "fdr"), pval]
  res[, log2OR:= log2(OR)]
  res <- res[, .(motif, cl, log2OR, padj)]
  if(plot)
  {
    vl_motif_cl_enrich_plot_only(obj = res,
                                 padj_cutoff= padj_cutoff,
                                 log2OR_cutoff= log2OR_cutoff,
                                 N_top= N_top,
                                 color_breaks= color_breaks,
                                 col= col,
                                 main= main,
                                 cex.balloons = cex.balloons, 
                                 auto_margins = auto_margins)
  }
  
  return(res)
}
