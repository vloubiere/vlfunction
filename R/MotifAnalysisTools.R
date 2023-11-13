#' Counts motifs
#'
#' Counts motif occurences in a set of regions
#'
#' @param sequences Named character vector of sequences to analyse. If provided takes over bed argument (in the case where both are specified)
#' @param bed Either a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param genome Genome to be used for coordinates ("dm6, "dm3") and as background for counting motifs when bg= "genome"
#' @param bg Background used to find motifs. Possible values include "genome" and "even". Default= "genome"
#' @param p.cutoff Pval cutoff used for motif detection
#' @param sel Vector of motif_ID to compute. see vl_Dmel_motifs_DB_full$motif_ID
#' @param motifDB The motifDB object to be used. see ?vl_Dmel_motifs_DB_full and ?vl_motifs_DB_v2, Default= vl_Dmel_motifs_DB_full
#' 
#' @examples
#' # Example run
#' sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", motif_ID]
#' 
#' # Regions
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
#' # Motif counts
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
                                      sel,
                                      genome,
                                      bg= "genome",
                                      p.cutoff= 5e-4,
                                      motifDB= vl_Dmel_motifs_DB_full)
{
  # Select motifs
  if(any(!sel %in% motifDB$motif_ID))
    stop("Some motif_ID(s) provided in 'sel' do not exist in motifDB$motif_ID")

  # Compute counts
  mot <- do.call(TFBSTools::PWMatrixList,
                 motifDB[match(unique(sel), motif_ID), pwms_log_odds])
  res <- as.matrix(motifmatchr::matchMotifs(mot,
                                            sequences,
                                            genome= genome,
                                            p.cutoff= p.cutoff,
                                            bg= bg,
                                            out= "scores")@assays@data[["motifCounts"]])
  
  # Save
  res <- as.data.table(res)
  setnames(res, as.character(sel))
  return(res)
}

#' Motif enrichment analysis
#' 
#' Compute motif enrichment between a set of regions and control regions
#'
#' @param counts data.table containing counts for the regions of interest
#' @param control.counts data.table containing counts for control regions (data.table)
#' @param names Convenient names to be used for plotting and so on... Default (NULL) returns to vl_Dmel_motifs_DB_full$motif_cluster names
#' @param plot Plot result?
#' @param padj.cutoff cutoff for plotting. Default to FALSE
#' @param top.enrich Show only n top enriched motifs
#' @param breaks Color breaks to be used. Defaults to range of filtered padj.
#' @param order Value to be used for ordering before selecting top enriched. Possible values are "padj", "log2OR". defaut= "padj"
#' @param add.motifs Should motif pwms be plotted?
#' @param cex.width expansion factor for motif widths
#' @param col Colors vector for bars
#' @param xlab x label
#' @param cex.height expansion factor for motif heights
#' 
#' @examples 
#' # Motif counts
#' sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", motif_ID]
#' top_SUHW <- vl_resizeBed(vl_SUHW_top_peaks,
#'                          upstream = 250,
#'                          downstream = 250,
#'                          genome = "dm3")
#' top_STARR <- vl_resizeBed(vl_STARR_DSCP_top_peaks,
#'                           upstream = 250,
#'                           downstream = 250,
#'                           genome = "dm3")
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
#' DT <- plot(DT,
#'            padj.cutoff= 1e-7)
#' vl_add_motifs(DT)
#' 
#' @return DT of enrichment values which can be plot using ?plot.vl_GO_enr
#' @export
vl_motif_enrich <- function(counts,
                            control.counts,
                            names= NULL,
                            plot= F,
                            padj.cutoff= 0.05,
                            top.enrich= NA, 
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
  
  # Append (cluster) name to motif ID ----
  if(is.null(names))
    names <- vl_Dmel_motifs_DB_full[names(obj)[-1], motif_cluster, on= "motif_ID"]
  names(obj)[-1] <- paste0(names, "__", names(obj)[-1]) # Append (cluster) name to motif ID
  
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
  res[, name:= tstrsplit(variable, "__", keep= 1)]
  res[, variable:= gsub(paste0("^", name, "__"), "", variable), name]
  
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
#' Compute motif enrichment for the cluster in cl_columns, using all the lines as background
#'
#' @param counts.list List of data.table containing counts for the regions of interest and potentially control regions (see next argument)
#' @param control.cl IDs of clusters to be used as background. default to NULL, meaning all clusters are used except the one being tested
#' @param names Convenient names to be used for plotting and so on... Default (NULL) returns vl_Dmel_motifs_DB_full$motif_cluster names
#' @param plot Should the result be plot using balloons plot? Default to FALSE
#' @param padj.cutoff cutoff for ballons to be plotted
#' @param top.enrich Select top enriched motifs/cluster. Default to NA (All)
#' @param order Value to be used for ordering before selecting top enriched. Possible values are "padj", "log2OR". defaut= "padj"
#' @param x.breaks Breaks used for ballon's sizes
#' @param color.breaks Color breaks used for coloring
#' @param col Vector of colors used for coloring
#' @param main Title. Default= NA
#' @param add.motifs Should motif pwms be plotted?
#' @param cex.width expansion factor for motif widths
#' @param cex.balloons Expansion factor for balloons
#' @param plot.empty.clusters Should empty clusters be plotted? Default= TRUE
#' @param cex.height expansion factor for motif heights
#'
#' @examples 
#' sel <- vl_Dmel_motifs_DB_full[collection=="jaspar", motif_ID]
#' top_SUHW <- vl_resizeBed(vl_SUHW_top_peaks, upstream = 250, downstream = 250, genome = "dm3")
#' top_STARR <- vl_resizeBed(vl_STARR_DSCP_top_peaks, upstream = 250, downstream = 250, genome = "dm3")
#' counts <- vl_motif_counts(top_SUHW, genome= "dm3", sel= sel)
#' control.counts <- vl_motif_counts(top_STARR, genome= "dm3", sel= sel)
#' counts.list <- list(SUHW= counts, ATAC= control.counts)
#' par(mar= c(3,20,2,5), las= 1)
#' DT <- vl_motif_cl_enrich(counts.list, padj.cutoff = 1e-3, plot= T, add.motifs= T)
#' DT <- plot(DT, padj.cutoff= 1e-5)
#' vl_add_motifs(DT)
#' 
#' @return Fisher test data.table.
#' @export
vl_motif_cl_enrich <- function(counts.list, 
                               control.cl= NULL,
                               names= NULL,
                               plot= F,
                               padj.cutoff= 1e-5,
                               top.enrich= NA,
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
#' @param sequences Named character vector of sequences to analyse. If provided takes over bed argument (in the case where both are specified)
#' @param bed Either a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param sel List of motifs to compute. see vl_Dmel_motifs_DB_full$motif
#' @param genome Genome to be used for coordinates ("dm6, "dm3") and as background for counting motifs when bg= "genome"
#' @param bg Background used to find motifs. Possible values include "genome" and "even". Default= "genome"
#' @param p.cutoff Pval cutoff used for motif detection
#' @param collapse.overlapping Should overlapping motifs be merged? If TRUE (default), motif instances that overlap more than 70 percent of their width are collapsed.
#' 
#' @examples
#' vl_motif_pos.data.table(vl_SUHW_top_peaks[1:2],
#'                         genome= "dm3",
#'                         sel= c("cisbp__M2328", "flyfactorsurvey__suHw_FlyReg_FBgn0003567", "jaspar__MA0533.1"))
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
                                   sel,
                                   genome,
                                   bg= "genome",
                                   p.cutoff= 5e-4,
                                   collapse.overlapping= TRUE,
                                   motifDB= vl_Dmel_motifs_DB_full)
{
  if(any(!sel %in% motifDB$motif_ID))
    stop("Some motif_ID(s) provided in 'sel' do not exist in motifDB$motif_ID")
  sub <- motifDB[match(unique(sel), motif_ID)]
  pos <- motifmatchr::matchMotifs(do.call(TFBSTools::PWMatrixList, 
                                          sub$pwms_log_odds),
                                  sequences,
                                  genome= genome,
                                  p.cutoff= p.cutoff,
                                  bg= bg,
                                  out= "positions")
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