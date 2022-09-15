#' Counts motifs
#'
#' Counts motif occurences in a set of regions
#'
#' @param sequences Named character vector of sequences to analyse. If provided takes over bed argument (in the case where both are specified)
#' @param bed Either a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param genome Genome to be used for coordinates ("dm6, "dm3")
#' @param p.cutoff Pval cutoff used for motif detection
#' @param sel List of motifs to compute. see vl_Dmel_motifs_DB_full$motif
#' @param what Waht value should be return. Either "motifMatches" (logical), "motifCounts" (numeric) or "motifScores" (numeric)
#' @examples
#' # Example run
#' suhw <- vl_motif_counts(vl_SUHW_top_peaks, "dm3")
#' ctls_regions <- vl_control_regions_BSgenome(bed= vl_SUHW_top_peaks, "dm3")
#' ctl <- vl_motif_counts(ctls_regions, "dm3")
#' pl <- vl_motif_enrich(suhw, ctl, plot= F)
#' par(mar= c(4,8,4,6))
#' plot(pl, padj_cutoff= 1e-3)
#' @return Matrix of motif counts
#' @export
vl_motif_counts <- function(sequences, ...) UseMethod("vl_motif_counts")

#' @describeIn vl_motif_counts Method to extract sequences from BSgenome
#' @export
vl_motif_counts.data.table <- function(bed, genome, ...)
{
  bed <- vl_importBed(bed)
  sequences <- BSgenome::getSeq(BSgenome::getBSgenome(genome), bed$seqnames, bed$start, bed$end, as.character= T)
  names(sequences) <- paste0(bed$seqnames, ":", bed$start, "-", bed$end, ":", bed$ranges)
  
  pl <- match.call()
  pl$bed <- pl$genome <- NULL
  pl$sequences <- sequences
  pl[[1]] <- quote(vl_motif_counts.character)
  eval(pl)
}

#' @describeIn vl_motif_counts Identify motifs in sequences
#' @export
vl_motif_counts.character <- function(sequences= NULL,
                                      p.cutoff= 5e-4,
                                      sel= vl_Dmel_motifs_DB_full[!is.na(vl_Dmel_motifs_DB_full$FBgn), motif],
                                      what= "motifCounts")
{
  if(is.null(names(sequences)))
    names(sequences) <- seq(sequences)
  
  # Select motifs
  if(any(!sel %in% vl_Dmel_motifs_DB_full$motif))
    stop("Some motif provided in 'sel' do not exist in vl_Dmel_motifs_DB_full$motif")
  sub <- vl_Dmel_motifs_DB_full[motif %in% sel]
  mot <- do.call(TFBSTools::PWMatrixList, 
                 sub$pwms_log_odds)
  
  # Count motifs
  res <- as.matrix(motifmatchr::matchMotifs(mot,
                                            sequences,
                                            p.cutoff= p.cutoff,
                                            bg= "even",
                                            out= "scores")@assays@data[[what]])
  colnames(res) <- sub$motif
  rownames(res) <- names(sequences)
  return(res)
}

#' Motif enrichment analysis
#' 
#' Compute motif enrichment between a set of regions and control regions
#'
#' @param counts matrix of counts for the regions of interest
#' @param control_counts matrix of counts for control regions
#' @param collapse_clusters A vector of clusters of lenght ncol(counts) used to collapse motifs (only the most enriched motif/cluster will be returned). Default= colnames(counts), meaning no collapsing
#' @param plot Plot result?
#' @param padj_cutoff cutoff for plotting
#' @param top_enrich Show only n top enriched motifs
#' @return DT of enrichment values
#' @examples 
#' # Example run
#' suhw <- vl_motif_counts(vl_SUHW_top_peaks, "dm3")
#' ctls_regions <- vl_control_regions_BSgenome(bed= vl_SUHW_top_peaks, "dm3")
#' ctl <- vl_motif_counts(ctls_regions, "dm3")
#' pl <- vl_motif_enrich(suhw, ctl, plot= F)
#' par(mar= c(4,8,4,6))
#' plot(pl, padj_cutoff= 1e-3)
#' @return Returns a table of enrichment which can be plot using ?plot.vl_GO_enr
#' @export
vl_motif_enrich <- function(counts,
                            control_counts,
                            collapse_clusters= colnames(counts),
                            plot= T,
                            padj_cutoff= 0.05,
                            top_enrich= Inf)
{
  if(!is.numeric(unlist(counts)))
    stop("counts_matrix should only contain numeric values")
  if(!is.numeric(unlist(control_counts)))
    stop("counts_matrix should only contain numeric values")
  if(length(collapse_clusters) != ncol(counts))
    stop("collapse_clusters shoulds match the lenght of ncol(counts)")
  
  # make obj
  obj <- rbindlist(list(set= as.data.table(counts),
                        control= as.data.table(control_counts)),
                   idcol = T)
  obj <- melt(obj, 
              id.vars = ".id",
              variable.name= "motif_ID")
  
  # Test enrichment
  res <- obj[, {
    # Contingency table 
    set <- factor(.id=="set", levels = c(T, F))
    motif <- factor(value>0, levels = c(T, F))
    tab <- table(set, motif)
    # Check contingency table -> Fisher
    res <- fisher.test(tab)
    .(OR= res$estimate,
      pval= res$p.value,
      set_hit= tab[1,1],
      set_total= sum(tab[1,]),
      ctl_hit= tab[2,1],
      ctl_total= sum(tab[2,]))
  }, motif_ID]
  
  # padj...
  res[, padj:= p.adjust(pval, method = "fdr"), pval]
  res[, log2OR:= log2(OR)]
  res$OR <- NULL
  
  # Collapsing
  res[, variable:= collapse_clusters[match(motif_ID, colnames(counts))], motif_ID]
  res <- res[, .SD[which.min(padj)], variable]
  
  # Order and save
  setcolorder(res,
              c("variable", "log2OR", "pval", "padj"))
  setattr(res, "class", c("vl_enr", "data.table", "data.frame"))
  
  if(plot)
    plot(res,
         padj_cutoff= padj_cutoff,
         top_enrich= top_enrich)
  
  invisible(res)
}

#' Compute motif enrichment
#'
#' Compute motif enrichment for the cluster in cl_columns, using all the lines as background
#'
#' @param counts_matrix A matrix containing motif counts (rows= sequences, cols= motifs)
#' @param cl_IDs vector of cluster IDs (used to split the data.table)
#' @param collapse_clusters A vector of clusters of lenght ncol(counts_matrix) used to collapse motifs (only the most enriched motif/cluster will be returned). Default= colnames(counts_matrix), meaning no collapsing
#' @param control_cl IDs of clusters to be used as background. default to unique(cl_IDs), meaning all clusters are used except the one being tested
#' @param plot Should the result be plot using balloons plot?
#' @param padj_cutoff cutoff for ballons to be plotted
#' @param top_enrich Select top enriched motifs/cluster
#' @param x_breaks Breaks used for ballon's sizes
#' @param color_breaks Color breaks used for coloring
#' @param col Vector of colors used for coloring
#' @param main Title. Default= NA
#' @examples 
#' # Sets of peaks
#' top_SUHW <- vl_resizeBed(vl_SUHW_top_peaks, upstream = 250, downstream = 250, genome = "dm3")
#' top_STARR <- vl_resizeBed(vl_STARR_DSCP_top_peaks, upstream = 250, downstream = 250, genome = "dm3")
#' rdm <- vl_random_regions_BSgenome("dm3", 2000, width= 500)
#' 
#' # Compute enrichment SUHW vs STARR-Seq
#' set <- rbind(top_SUHW,
#'              top_STARR, fill= T)
#' counts <- vl_motif_counts(set, "dm3")
#' enr <- vl_motif_cl_enrich(counts, 
#'                           cl_IDs = c(rep(1, nrow(top_SUHW)),
#'                                      rep(2, nrow(top_STARR))),
#'                           plot=F)
#' par(las= 1)
#' plot(enr, padj_cutoff= 1e-5)
#' 
#' Compute enrichment of SUHW & STARR-Seq over random regions
#' counts <- vl_motif_counts(rbind(top_SUHW,
#'                                 top_STARR, 
#'                                 rdm, fill= T), "dm3")
#' cl_IDs <- c(rep("suhw", nrow(top_SUHW)),
#'             rep("STARR", nrow(top_STARR)),
#'             rep("rdm", nrow(rdm)))
#' enr <- vl_motif_cl_enrich(counts, 
#'                           cl_IDs = cl_IDs,
#'                           control_cl = "rdm",
#'                           plot=F)
#' par(las= 1)
#' plot(enr, padj_cutoff= 1e-20, top_enrich= 20)
#' @return Fisher test data.table.
#' @export
vl_motif_cl_enrich <- function(counts_matrix, 
                               cl_IDs,
                               collapse_clusters= colnames(counts_matrix),
                               control_cl= unique(cl_IDs),
                               plot= T,
                               padj_cutoff= 1e-5,
                               top_enrich= Inf,
                               x_breaks,
                               color_breaks,
                               cex.balloons= 1,
                               col= c("cornflowerblue", "lightgrey", "tomato"),
                               main= NA)
{
  if(!is.numeric(unlist(counts_matrix)))
    stop("counts_matrix should only contain numeric values")
  if(!is.data.table(counts_matrix))
    counts_matrix <- as.data.table(counts_matrix)
  if(!is.factor(cl_IDs))
    cl_IDs <- factor(cl_IDs)

  # Compute enrichment
  res <- lapply(levels(cl_IDs), function(cl) 
  {
    vl_motif_enrich(counts = counts_matrix[cl_IDs==cl],
                    control_counts = counts_matrix[cl_IDs!=cl & cl_IDs %in% control_cl],
                    collapse_clusters= collapse_clusters,
                    plot= F)
  })
  names(res) <- levels(cl_IDs)
  res <- rbindlist(res, idcol = "cl")
  res[, cl:= factor(cl, unique(cl))]
  class(res) <- c("vl_enr_cl", "data.table", "data.frame")
  
  # plot
  if(plot)
  {
    pl <- match.call()
    pl$counts_matrix <- pl$cl_IDs <- pl$control_cl <- pl$plot <- NULL
    pl$obj <- res
    pl[[1]] <- quote(plot)
    eval(pl)
  }
  
  return(res)
}