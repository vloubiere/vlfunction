#' Counts motifs
#'
#' Counts motif occurences in a set of regions
#'
#' @param sequences Named character vector of sequences to analyse. If provided takes over bed argument (in the case where both are specified)
#' @param bed Either a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start' and 'end' columns.
#' @param sel Vector of motif_ID(s) existing in motifDB (see below) for which motif matches should be counted. Defaults to jaspar PWMs.
#' @param motifDB The motifDB to be used (see 'sel' argument; ?vl_Dmel_motifs_DB_full and ?vl_motifs_DB_v2 for details). Default= vl_Dmel_motifs_DB_full.
#' @param pwms_log_odds A PWMatrixList (in log2 odd ratio format) for which motif matches should be counted. Overrides sel and motifDB arguments (see above).
#' @param genome Genome to be used for coordinates ("dm6, "dm3") and as background for counting motifs when bg= "genome".
#' @param bg Background used to find motifs. Possible values include "genome" and "even". Default= "genome"
#' @param p.cutoff p.value cutoff used for motif detection. Default= 5e-4.
#' 
#' @examples
#' # Select motifs
#' sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", motif_ID]
#' 
#' # Resize example peaks
#' top_SUHW <- vl_resizeBed(vl_SUHW_top_peaks,
#'                          upstream = 250,
#'                          downstream = 250,
#'                          genome = "dm3")
#' top_STARR <- vl_resizeBed(vl_STARR_DSCP_top_peaks,
#'                           upstream = 250,
#'                           downstream = 250,
#'                           genome = "dm3")
#' ctls_regions <- vl_control_regions_BSgenome(bed= vl_SUHW_top_peaks,
#'                                             genome= "dm3")
#'                           
#' # Count motifs
#' suhw <- vl_motif_counts(top_SUHW,
#'                         genome= "dm3",
#'                         sel= sel)
#' starr <- vl_motif_counts(top_STARR,
#'                          genome= "dm3",
#'                          sel= sel)
#' ctl <- vl_motif_counts(ctls_regions,
#'                        genome= "dm3",
#'                        sel= sel)
#'                        
#' # Enrichment                        
#' pl <- vl_motif_enrich(suhw,
#'                       ctl,
#'                       plot= F)
#'
#' # Plot
#' plot(pl,
#'      padj.cutoff= 5e-2)
#' vl_motif_enrich(starr, ctl)  
#' 
#' @return Matrix of motif counts
#' @export
vl_motif_counts <- function(sequences, ...) UseMethod("vl_motif_counts")

#' @describeIn vl_motif_counts Method to extract sequences from BSgenome
#' @export
vl_motif_counts.data.table <- function(bed,
                                       genome,
                                       ...)
{
  sequences <- vl_getSequence(bed, genome)
  vl_motif_counts.character(sequences,
                            genome= genome,
                            ...)
}

#' @describeIn vl_motif_counts Identify motifs in sequences
#' @export
vl_motif_counts.character <- function(sequences= NULL,
                                      sel= vl_Dmel_motifs_DB_full[collection=="jaspar", motif_ID],
                                      motifDB= vl_Dmel_motifs_DB_full,
                                      pwm_log_odds= NULL,
                                      genome,
                                      bg= "genome",
                                      p.cutoff= 5e-4)
{
  # Select motifs
  if(is.null(pwm_log_odds) && any(!sel %in% motifDB$motif_ID))
    stop("Some motif_ID(s) provided in 'sel' do not exist in motifDB$motif_ID.")
  if(is.null(pwm_log_odds) && length(sel)!=length(unique(sel)))
    stop("Selected motif_ID(s) should be unique.")
  if(is.null(pwm_log_odds))
    pwm_log_odds <- do.call(TFBSTools::PWMatrixList,
                            motifDB[match(sel, motif_ID), pwms_log_odds])

  # Compute counts
  res <- motifmatchr::matchMotifs(pwm_log_odds,
                                  sequences,
                                  genome= genome,
                                  p.cutoff= p.cutoff,
                                  bg= bg,
                                  out= "scores")@assays@data[["motifCounts"]]
  res <- as.matrix(res)
  res <- as.data.table(res)
  setnames(res,
           TFBSTools::name(pwm_log_odds))
  
  # Save
  return(res)
}

#' Motif enrichment analysis
#' 
#' Compute motif enrichment between a set of regions and control regions
#'
#' @param counts data.table containing counts for the regions of interest.
#' @param control.counts data.table containing counts for control regions. Should have the same number of columns as 'counts' table.
#' @param names Convenient names to be used for aggregating and plotting. By default, returns motif_cluster matching motif_ID in vl_Dmel_motifs_DB_full.
#' @param plot Plot result?
#' @param padj.cutoff cutoff for plotting. Default= FALSE.
#' @param top.enrich Show only n top enriched motifs. Default= Inf (all).
#' @param min.counts The minimum number of counts required to call a hit. Default= 3L.
#' @param breaks Color breaks to be used. Defaults to range of filtered padj.
#' @param order Value to be used for ordering before selecting top enriched. Possible values are "padj", "log2OR". Defaut= "padj".
#' @param add.motifs Should motif pwms be plotted?
#' @param cex.width Expansion factor for motif widths.
#' @param col Colors vector for bars.
#' @param xlab x label.
#' @param cex.height Expansion factor for motif heights.
#' 
#' @examples 
#' # Select motifs
#' sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", motif_ID]
#' 
#' # Resize example peaks
#' top_SUHW <- vl_resizeBed(vl_SUHW_top_peaks,
#'                          upstream = 250,
#'                          downstream = 250,
#'                          genome = "dm3")
#' top_STARR <- vl_resizeBed(vl_STARR_DSCP_top_peaks,
#'                           upstream = 250,
#'                           downstream = 250,
#'                           genome = "dm3")
#'                           
#' # Count motifs
#' counts <- vl_motif_counts(top_SUHW,
#'                           genome= "dm3",
#'                           sel= sel)
#' control.counts <- vl_motif_counts(top_STARR,
#'                                   genome= "dm3",
#'                                   sel= sel)
#'                                   
#' # Enrichment
#' DT <- vl_motif_enrich(counts,
#'                       control.counts,
#'                       add.motifs= T,
#'                       padj.cutoff= 1e-5,
#'                       plot= T)
#' 
#' # Plot
#' par(mai= c(.9, 2, .9, .9))
#' DT <- plot(DT,
#'            padj.cutoff= 1e-7)
#' vl_add_motifs(DT)
#' 
#' @return DT of enrichment values which can be plot using ?plot.vl_GO_enr
#' @export
vl_motif_enrich <- function(counts,
                            control.counts,
                            names= vl_Dmel_motifs_DB_full[colnames(counts), motif_cluster, on= "motif_ID"],
                            plot= F,
                            padj.cutoff= 0.05,
                            top.enrich= Inf,
                            min.counts= 3L,
                            order= "padj",
                            breaks= NULL,
                            col= c("blue", "red"),
                            xlab = "Odd Ratio (log2)",
                            add.motifs= F,
                            cex.width= 1,
                            cex.height= 1)
{
  if(!is.data.table(counts))
    counts <- as.data.table(counts)
  if(!is.data.table(control.counts))
    counts <- as.data.table(control.counts)
  if(!all(sapply(counts, class)=="numeric"))
    stop("counts should only contain numeric values")
  if(!all(sapply(control.counts, class)=="numeric"))
    stop("control.counts should only contain numeric values")
  if(!is.null(names) && length(names)!=ncol(counts))
    stop("names should match ncol(counts)")
  
  # make obj ----
  obj <- rbindlist(list(set= counts,
                        control= control.counts),
                   idcol = T)
  
  # Melt ----
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
  res[, name:= names]
  res[is.na(name), name:= variable] # If names can't be found
  
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
                      padj.cutoff= padj.cutoff,
                      top.enrich= top.enrich,
                      min.counts= min.counts,
                      order= order,
                      breaks= breaks, 
                      xlab = xlab, 
                      col = col)
    if(add.motifs)
      vl_add_motifs(DT,
                    cex.width= cex.width, 
                    cex.height= cex.height)
  }
  invisible(res)
}

#' Compute motif enrichment
#'
#' Compute motif enrichment for the cluster in cl_columns, using all the lines as background.
#'
#' @param counts.list List of data.table containing counts for the regions of interest and potentially control regions (see next argument).
#' @param control.cl IDs of clusters to be used as background. Default= NULL, meaning all clusters are used except the one being tested.
#' @param names Convenient names to be used for aggregating and plotting. By default, returns motif_cluster matching motif_ID in vl_Dmel_motifs_DB_full.
#' @param plot Should the result be plotted using balloons plot? Default to FALSE
#' @param padj.cutoff Cutoff for balloons to be plotted.
#' @param top.enrich Select top enriched motifs/cluster. Default= Inf (All).
#' @param min.counts The minimum number of counts required to call a hit. Default= 3L.
#' @param order Value to be used for ordering before selecting top enriched. Possible values are "padj" (the default) and "log2OR".
#' @param x.breaks Breaks used for balloon's sizes.
#' @param color.breaks Color breaks used for coloring.
#' @param col Vector of colors used for coloring.
#' @param main Title. Default= NA.
#' @param add.motifs Should motif pwms be plotted?
#' @param cex.width Expansion factor for motif widths.
#' @param cex.balloons Expansion factor for balloons.
#' @param plot.empty.clusters Should empty clusters be plotted? Default= TRUE.
#' @param cex.height Expansion factor for motif heights.
#'
#' @examples 
#' # Select motifs
#' sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", motif_ID]
#' 
#' # Resize example peaks
#' top_SUHW <- vl_resizeBed(vl_SUHW_top_peaks, upstream = 250, downstream = 250, genome = "dm3")
#' top_STARR <- vl_resizeBed(vl_STARR_DSCP_top_peaks, upstream = 250, downstream = 250, genome = "dm3")
#' 
#' # Count motifs
#' counts <- vl_motif_counts(top_SUHW, genome= "dm3", sel= sel)
#' control.counts <- vl_motif_counts(top_STARR, genome= "dm3", sel= sel)
#' counts.list <- list(SUHW= counts,
#'                     ATAC= control.counts)
#' 
#' # Enrichment and plot
#' vl_par(mai= c(.9, 2, .9, .9))
#' DT <- vl_motif_cl_enrich(counts.list, padj.cutoff = 1e-3, plot= T, add.motifs= T)
#' DT <- plot(DT, padj.cutoff= 1e-5)
#' vl_add_motifs(DT)
#' 
#' @return Fisher test data.table.
#' @export
vl_motif_cl_enrich <- function(counts.list, 
                               control.cl= NULL,
                               names= vl_Dmel_motifs_DB_full[colnames(counts), motif_cluster, on= "motif_ID"],
                               plot= F,
                               padj.cutoff= 1e-5,
                               top.enrich= Inf,
                               min.counts= 3L,
                               order= "padj",
                               x.breaks,
                               color.breaks,
                               cex.balloons= 1,
                               col= c("cornflowerblue", "lightgrey", "tomato"),
                               main= NA,
                               plot.empty.clusters= T,
                               add.motifs= F,
                               cex.width= 1,
                               cex.height= 1)
{
  if(!all(sapply(counts.list, is.data.table)))
    counts.list <- lapply(counts.list, as.data.table)
  if(is.null(names(counts.list)))
    names(counts.list) <- seq(counts.list)
  if(!is.null(control.cl) && any(!control.cl %in% names(counts.list)))
    stop("control.cl should match names(counts.list)")
  
  # Compute enrichment in non-control clusters
  cmb <- data.table(cl= names(counts.list))
  if(is.null(control.cl))
    cmb <- cmb[, .(ctl= cmb$cl[cmb$cl!=cl]), cl] else
      cmb <- cmb[, .(ctl= control.cl), cl]
  # Remove self-comparisons ----
  cmb[, self:= identical(cl, ctl), cl]
  cmb <- cmb[!(self), !"self"]
  # Motif enrichment ----
  res <- cmb[, {
    vl_motif_enrich(counts = counts.list[[cl]],
                    control.counts = rbindlist(counts.list[ctl]),
                    names= names,
                    plot= F)
  }, cl]
  res[, cl:= factor(cl, unique(cl))]
  setattr(res, "class", c("vl_enr_cl", "data.table", "data.frame"))

  # plot
  if(plot)
  {
    DT <- plot.vl_enr_cl(obj = res,
                         padj.cutoff= padj.cutoff,
                         top.enrich= top.enrich, 
                         min.counts= min.counts,
                         order= order,
                         x.breaks= x.breaks,
                         color.breaks= color.breaks,
                         cex.balloons= cex.balloons,
                         col= col,
                         main= main, 
                         plot.empty.clusters = plot.empty.clusters)
    if(add.motifs)
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
#' @param sequences Named character vector of sequences to analyse. If provided takes over bed argument (in the case where both are specified).
#' @param bed Either a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start' and 'end' columns.
#' @param sel Vector of motif_ID(s) existing in motifDB (see below) and for which motif matches should be mapped Defaults to jaspar PWMs.
#' @param motifDB The motifDB to be used (see 'sel' argument; ?vl_Dmel_motifs_DB_full and ?vl_motifs_DB_v2 for details). Default= vl_Dmel_motifs_DB_full.
#' @param pwms_log_odds A PWMatrixList (in log2 odd ratio format) for which motif matches should be mapped. Overrides sel and motifDB arguments (see above).
#' @param genome Genome to be used for coordinates ("dm6, "dm3") and as background for counting motifs when bg= "genome".
#' @param bg Background used to find motifs. Possible values include "genome" and "even". Default= "genome"
#' @param p.cutoff p.value cutoff used for motif detection. Default= 5e-4.
#' @param collapse.overlapping Should overlapping motifs be merged? If TRUE (default), motif instances that overlap more than 70 percent of their width are collapsed.
#' 
#' @examples
#' test <- vl_motif_pos.data.table(vl_SUHW_top_peaks[1:2],
#'                                 genome= "dm3",
#'                                 sel= c("cisbp__M2328",
#'                                        "flyfactorsurvey__suHw_FlyReg_FBgn0003567",
#'                                        "jaspar__MA0533.1"))
#'                         
#' @return A list of positions of length = length(sequences) 
#' @export
vl_motif_pos <- function(sequences, ...) UseMethod("vl_motif_pos")

#' @describeIn vl_motif_pos Compute motifs position within
#' @export
vl_motif_pos.data.table <- function(bed, genome, ...)
{
  sequences <- vl_getSequence(bed, genome)
  vl_motif_pos.character(sequences, genome= genome, ...)
}

#' @describeIn vl_motif_pos Identify motif positions within sequences
#' @export
vl_motif_pos.character <- function(sequences,
                                   sel= vl_Dmel_motifs_DB_full[collection=="jaspar", motif_ID],
                                   motifDB= vl_Dmel_motifs_DB_full,
                                   pwm_log_odds= NULL,
                                   genome,
                                   bg= "genome",
                                   p.cutoff= 5e-4,
                                   collapse.overlapping= TRUE)
{
  # Checks
  if(is.null(pwm_log_odds) && any(!sel %in% motifDB$motif_ID))
    stop("Some motif_ID(s) provided in 'sel' do not exist in motifDB$motif_ID.")
  if(is.null(pwm_log_odds) && length(sel)!=length(unique(sel)))
    stop("Selected motif_ID(s) should be unique.")
  if(is.null(pwm_log_odds))
    pwm_log_odds <- do.call(TFBSTools::PWMatrixList,
                            motifDB[match(sel, motif_ID), pwms_log_odds])
  
  # Map motifs
  pos <- motifmatchr::matchMotifs(pwm_log_odds,
                                  sequences,
                                  genome= genome,
                                  p.cutoff= p.cutoff,
                                  bg= bg,
                                  out= "positions")
  # Name
  names(pos) <- TFBSTools::name(pwm_log_odds)
  
  # Format and collapse if specified
  pos <- lapply(pos, function(x)
  {
    names(x) <- seq(length(sequences))
    lapply(x, function(y) {
      y <- as.data.table(y)
      # Collapse motifs that overlap more than 70% -> return merged motifs
      if(collapse.overlapping && nrow(y))
      {
        # y <- data.table(start=c(1,2,3,4,5), end= c(5,6,7,8,9), score= 1, width= 4)
        setorderv(y, c("start", "end"))
        motif_width <- ceiling(y[1, width]*0.7)
        # Identify overlapping motifs
        y[, idx:= cumsum(c(1, (start[-.N]+motif_width-1)<=start[-1]))]
        # Re-split tightly clustered motifs 
        y[, idx2:= round((start-start[1])/motif_width[1]), idx]
        y$idx <- rleidv(y, c("idx", "idx2"))
        # Return motifs
        y <- y[, .(start= start[1], 
                   end= end[.N],
                   score= max(score)), idx]
        y[, width:= end-start+1]
        y[, .(start, end, width, score)]
      }
      return(y)
      })
  })
  return(pos)
}

#' PWM perc to log2
#'
#' @param perc_pwm A percentage PWM, where all columns sum to 1
#'
#' @return Returns a log2 odd ratio pwm
#' @export
#'
#' @examples
#' vl_Dmel_motifs_DB_full$pwms_perc[[5]]@profileMatrix
#' vl_pwm_perc_to_log2(vl_Dmel_motifs_DB_full$pwms_perc[[5]]@profileMatrix)
#' 
vl_pwm_perc_to_log2 <- function(perc_pwm)
{
  log2(perc_pwm/matrix(0.25, nrow= 4, ncol= ncol(perc_pwm))+9.9999e-06)
}
