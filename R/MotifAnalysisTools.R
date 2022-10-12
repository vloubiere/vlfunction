#' Counts motifs
#'
#' Counts motif occurences in a set of regions
#'
#' @param sequences Named character vector of sequences to analyse. If provided takes over bed argument (in the case where both are specified)
#' @param bed Either a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param genome Genome to be used for coordinates ("dm6, "dm3")
#' @param p.cutoff Pval cutoff used for motif detection
#' @param sel Vector of motif_ID to compute. see vl_Dmel_motifs_DB_full$motif_ID
#' @param what Waht value should be return. Either "motifMatches" (logical), "motifCounts" (numeric) or "motifScores" (numeric)
#' @examples
#' # Example run
#' sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", motif_ID]
#' top_SUHW <- vl_resizeBed(vl_SUHW_top_peaks, upstream = 250, downstream = 250, genome = "dm3")
#' top_STARR <- vl_resizeBed(vl_STARR_DSCP_top_peaks, upstream = 250, downstream = 250, genome = "dm3")
#' counts <- vl_motif_counts(top_SUHW, genome= "dm3", sel= sel)
#' control_counts <- vl_motif_counts(top_STARR, genome= "dm3", sel= sel)
#' 
#' suhw <- vl_motif_counts(vl_SUHW_top_peaks, "dm3")
#' ctls_regions <- vl_control_regions_BSgenome(bed= vl_SUHW_top_peaks, "dm3")
#' ctl <- vl_motif_counts(ctls_regions, "dm3")
#' pl <- vl_motif_enrich(suhw, ctl, plot= F)
#' par(mar= c(4,15,4,6))
#' plot(pl, padj_cutoff= 1e-3)
#' @return Matrix of motif counts
#' @export
vl_motif_counts <- function(sequences, ...) UseMethod("vl_motif_counts")

#' @describeIn vl_motif_counts Method to extract sequences from BSgenome
#' @export
vl_motif_counts.data.table <- function(bed, genome, ...)
{
  sequences <- vl_getSequence(bed, genome)
  vl_motif_counts.character(sequences, ...)
}

#' @describeIn vl_motif_counts Identify motifs in sequences
#' @export
vl_motif_counts.character <- function(sequences= NULL,
                                      p.cutoff= 5e-4,
                                      sel,
                                      what= "motifCounts")
{
  # Select motifs
  if(any(!sel %in% vl_Dmel_motifs_DB_full$motif_ID))
    stop("Some motif_ID(s) provided in 'sel' do not exist in vl_Dmel_motifs_DB_full$motif_ID")
  mot <- do.call(TFBSTools::PWMatrixList, 
                 vl_Dmel_motifs_DB_full[match(unique(sel), motif_ID), pwms_log_odds])
  # Count
  res <- as.matrix(motifmatchr::matchMotifs(mot,
                                            sequences,
                                            p.cutoff= p.cutoff,
                                            bg= "even",
                                            out= "scores")@assays@data[[what]])
  res <- as.data.table(res)
  setnames(res, as.character(sel))
  return(res)
}

#' Motif enrichment analysis
#' 
#' Compute motif enrichment between a set of regions and control regions
#'
#' @param counts data.table containing counts for the regions of interest
#' @param control_counts data.table containing counts for control regions (data.table)
#' @param names Convenient names to be used for plotting and so on... Default (NULL) returns to vl_Dmel_motifs_DB_full$motif_cluster
#' @param plot Plot result?
#' @param padj_cutoff cutoff for plotting. Default to FALSE
#' @param top_enrich Show only n top enriched motifs
#' @param breaks Color breaks to be used. Defaults to range of filtered padj.
#' @param order Value to be used for ordering before selecting top enriched. Possible values are "padj", "log2OR". defaut= "padj"
#' @param add_motifs Should motif pwms be plotted?
#' @param cex.width expansion factor for motif widths
#' @param col Colors vector for bars
#' @param xlab x label
#' @param cex.height expansion factor for motif heights
#'
#' @return DT of enrichment values
#' @examples 
#' # Example run
#' sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", motif_ID]
#' top_SUHW <- vl_resizeBed(vl_SUHW_top_peaks, upstream = 250, downstream = 250, genome = "dm3")
#' top_STARR <- vl_resizeBed(vl_STARR_DSCP_top_peaks, upstream = 250, downstream = 250, genome = "dm3")
#' counts <- vl_motif_counts(top_SUHW, genome= "dm3", sel= sel)
#' control_counts <- vl_motif_counts(top_STARR, genome= "dm3", sel= sel)
#' par(mar= c(4,15,4,6))
#' DT <- vl_motif_enrich(counts, control_counts, add_motifs= T, padj_cutoff= 1e-5, plot= T)
#' DT <- plot(DT, padj_cutoff= 1e-7)
#' vl_add_motifs(DT)
#' @return Returns a table of enrichment which can be plot using ?plot.vl_GO_enr
#' @export
vl_motif_enrich <- function(counts,
                            control_counts,
                            names= NULL,
                            plot= F,
                            padj_cutoff= 0.05,
                            top_enrich= NA, 
                            order= "padj",
                            breaks= NULL,
                            col= c("blue", "red"),
                            xlab = "Odd Ratio (log2)",
                            add_motifs= F,
                            cex.width= 1,
                            cex.height= 1)
{
  if(!is.data.table(counts))
    counts <- as.data.table(counts)
  if(!is.data.table(control_counts))
    counts <- as.data.table(control_counts)
  if(!all(sapply(counts, class)=="numeric"))
    stop("counts should only contain numeric values")
  if(!all(sapply(control_counts, class)=="numeric"))
    stop("control_counts should only contain numeric values")
  if(!is.null(names) && length(names)!=ncol(counts))
    stop("names should match ncol(counts)")
  
  # make obj
  obj <- rbindlist(list(set= as.data.table(counts),
                        control= as.data.table(control_counts)),
                   idcol = T)
  obj <- melt(obj, 
              id.vars = ".id",
              variable.name= "variable")
  
  # Test enrichment
  res <- obj[, {
    # Contingency table
    tab <- table(factor(.id=="set", levels = c(T, F)), 
                 factor(value>0, levels = c(T, F)))
    # Check contingency table -> Fisher
    res <- fisher.test(tab)
    .(OR= res$estimate,
      pval= res$p.value,
      set_hit= sum(.id=="set" & value>0),
      set_total= sum(.id=="set"),
      ctl_hit= sum(.id=="control" & value>0),
      ctl_total= sum(.id=="control"))
  }, variable]
  
  # Add names
  if(is.null(names))
    res[vl_Dmel_motifs_DB_full, name:= motif_cluster, on= "variable==motif_ID"] else
      res[, names:= names]
  
  # padj...
  res[, padj:= p.adjust(pval, method = "fdr")]
  res[, log2OR:= log2(OR)]
  res$OR <- NULL
  
  # Order and save
  setcolorder(res, c("variable", "log2OR", "pval", "padj"))
  setattr(res, "class", c("vl_enr", "data.table", "data.frame"))
  if(plot)
  {
    DT <- plot.vl_enr(obj= res,
                      padj_cutoff= padj_cutoff,
                      top_enrich= top_enrich,
                      order= order,
                      breaks= breaks, 
                      xlab = xlab, 
                      col = col)
    if(add_motifs)
      vl_add_motifs(DT,
                    cex.width= cex.width, 
                    cex.height= cex.height)
  }
  invisible(res)
}

#' Compute motif enrichment
#'
#' Compute motif enrichment for the cluster in cl_columns, using all the lines as background
#'
#' @param counts_list List of data.table containing counts for the regions of interest and potentially control regions
#' @param control_cl IDs of clusters to be used as background. default to NULL, meaning all clusters are used except the one being tested
#' @param names Convenient names to be used for plotting and so on... Default (NULL) returns vl_Dmel_motifs_DB_full$motif_cluster
#' @param plot Should the result be plot using balloons plot? Default to FALSE
#' @param padj_cutoff cutoff for ballons to be plotted
#' @param top_enrich Select top enriched motifs/cluster. Default to NA (All)
#' @param order Value to be used for ordering before selecting top enriched. Possible values are "padj", "log2OR". defaut= "padj"
#' @param x_breaks Breaks used for ballon's sizes
#' @param color_breaks Color breaks used for coloring
#' @param col Vector of colors used for coloring
#' @param main Title. Default= NA
#' @param add_motifs Should motif pwms be plotted?
#' @param cex.width expansion factor for motif widths
#' @param cex.balloons Expansion factor for balloons
#' @param plot_empty_clusters Should empty clusters be plotted? Default= TRUE
#' @param cex.height expansion factor for motif heights
#'
#' @examples 
#' sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", variable]
#' top_SUHW <- vl_resizeBed(vl_SUHW_top_peaks, upstream = 250, downstream = 250, genome = "dm3")
#' top_STARR <- vl_resizeBed(vl_STARR_DSCP_top_peaks, upstream = 250, downstream = 250, genome = "dm3")
#' counts <- vl_motif_counts(top_SUHW, genome= "dm3", sel= sel)
#' control_counts <- vl_motif_counts(top_STARR, genome= "dm3", sel= sel)
#' counts_list <- list(SUHW= counts, ATAC= control_counts)
#' par(mar= c(3,20,2,5), las= 1)
#' DT <- vl_motif_cl_enrich(counts_list, padj_cutoff = 1e-3, plot= T, add_motifs= T)
#' DT <- plot(DT, padj_cutoff= 1e-5)
#' vl_add_motifs(DT)
#' 
#' @return Fisher test data.table.
#' @export
vl_motif_cl_enrich <- function(counts_list, 
                               control_cl= NULL,
                               names= NULL,
                               plot= F,
                               padj_cutoff= 1e-5,
                               top_enrich= NA,
                               order= "padj",
                               x_breaks,
                               color_breaks,
                               cex.balloons= 1,
                               col= c("cornflowerblue", "lightgrey", "tomato"),
                               main= NA,
                               plot_empty_clusters= T,
                               add_motifs= F,
                               cex.width= 1,
                               cex.height= 1)
{
  if(!all(sapply(counts_list, is.data.table)))
    counts_list <- lapply(counts_list, as.data.table)
  if(is.null(names(counts_list)))
    names(counts_list) <- seq(counts_list)
  if(!is.null(control_cl) && any(!control_cl %in% names(counts_list)))
    stop("control_cl should match names(counts_list)")
  
  # Compute enrichment in non-control clusters
  cmb <- CJ(cl= names(counts_list), 
            ctl= names(counts_list), 
            unique = T, 
            sorted = F)
  cmb <- cmb[cl!=ctl]
  if(!is.null(control_cl))
    cmb <- cmb[!(cl %in% control_cl) & (ctl %in% control_cl)]
  res <- cmb[, {
    vl_motif_enrich(counts = counts_list[[cl]],
                    control_counts = rbindlist(counts_list[ctl]),
                    names= names,
                    plot= F)
  }, cl]
  res[, cl:= factor(cl, unique(cl))]
  setattr(res, "class", c("vl_enr_cl", "data.table", "data.frame"))

  # plot
  if(plot)
  {
    DT <- plot.vl_enr_cl(obj = res,
                         padj_cutoff= padj_cutoff,
                         top_enrich= top_enrich, 
                         order= order,
                         x_breaks= x_breaks,
                         color_breaks= color_breaks,
                         cex.balloons= cex.balloons,
                         col= col,
                         main= main, 
                         plot_empty_clusters = plot_empty_clusters)
    if(add_motifs)
      vl_add_motifs(DT,
                    cex.width= cex.width, 
                    cex.height= cex.height)
  }
  
  return(res)
}

#' Find motif positions
#'
#' Compute motif positions
#'
#' @param sequences Named character vector of sequences to analyse. If provided takes over bed argument (in the case where both are specified)
#' @param bed Either a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param genome Genome to be used for coordinates ("dm6, "dm3")
#' @param p.cutoff Pval cutoff used for motif detection
#' @param sel List of motifs to compute. see vl_Dmel_motifs_DB_full$motif
#' @examples
#' vl_motif_pos.data.table(vl_SUHW_top_peaks[1:2], "dm3", sel= c("cisbp__M2328", "flyfactorsurvey__suHw_FlyReg_FBgn0003567", "jaspar__MA0533.1"))
#' @return A list of positions of lenght = length(sequences) 
#' @export
vl_motif_pos <- function(sequences, ...) UseMethod("vl_motif_pos")

#' @describeIn vl_motif_pos Compute motifs position within
#' @export
vl_motif_pos.data.table <- function(bed, genome, ...)
{
  sequences <- vl_getSequence(bed, genome)
  vl_motif_pos.character(sequences, ...)
}

#' @describeIn vl_motif_pos Identify motif positions within sequences
#' @export
vl_motif_pos.character <- function(sequences,
                                   p.cutoff= 5e-4,
                                   sel)
{
  if(any(!sel %in% vl_Dmel_motifs_DB_full$motif_ID))
    stop("Some motif_ID(s) provided in 'sel' do not exist in vl_Dmel_motifs_DB_full$motif_ID")
  sub <- vl_Dmel_motifs_DB_full[match(unique(sel), motif_ID)]
  pos <- motifmatchr::matchMotifs(do.call(TFBSTools::PWMatrixList, 
                                          sub$pwms_log_odds),
                                  sequences,
                                  p.cutoff= p.cutoff,
                                  bg= "even",
                                  out= "positions")
  pos <- lapply(pos, function(x)
  {
    names(x) <- seq(length(sequences))
    lapply(x, function(y) as.data.table(y))
  })
  return(pos)
}


