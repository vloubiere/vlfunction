#' Counts motifs
#'
#' Counts motif occurences in a set of regions
#'
#' @param sequences Named character vector of sequences to analyse. If provided takes over bed argument (in the case where both are specified)
#' @param bed Either a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start' and 'end' columns.
#' @param pwms_log_odds A PWMatrixList (in log2 odd ratio format) for which motif matches should be counted. Overrides sel and motifDB arguments (see above).
#' @param genome Genome to be used for coordinates ("dm6, "dm3") and as background for counting motifs when bg= "genome".
#' @param bg Background used to find motifs. Possible values include "genome" and "even". Default= "genome"
#' @param p.cutoff p.value cutoff used for motif detection. For enrichment analyses based on presence/absence of a motif, high cutoff might perform better (1e-4 or 5e-5) while for regression analyses, lower cutoffs might be preferred (5e-4). Default= 5e-5 (stringent).
#' 
#' @examples
#' # Resize example peaks
#' SUHW <- vl_resizeBed(vl_SUHW_top_peaks, genome = "dm3")
#' STARR <- vl_resizeBed(vl_STARR_DSCP_top_peaks, genome = "dm3")
#' 
#' # Generate same number of random regions
#' random <- vl_control_regions_BSgenome(bed= STARR, genome= "dm3")
#' 
#' # Count JAPSPAR motifs (see below to use custom list of PWMs)
#' suhw <- vl_motif_counts(SUHW, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' starr <- vl_motif_counts(top_STARR, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' ctl <- vl_motif_counts(random, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' 
#' # Starting from sequence instead of bed file
#' seq <- vl_getSequence(SUHW, genome= "dm3")
#' seq_suhw <- vl_motif_counts(seq, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' identical(suhw, seq_suhw)
#' 
#' # Motifs can also be counted using a custom PWMatrixList, for example for promoter motifs:
#' prom_db <- readRDS("/groups/stark/almeida/data/motifs/CP_motifs/CP_motifs_PWM.rds")
#' prom <- prom_db$Pwms_log_odds
#' for(i in seq(prom))
#'   prom[[i]]@profileMatrix <- vl_pwm_perc_to_log2(prom_db$Pwms_perc[[i]]@profileMatrix)
#' 
#' prom_motifs <- vl_motif_counts(STARR,
#'                                pwm_log_odds= prom,
#'                                genome= "dm3")
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
  vl_motif_counts.default(sequences,
                          genome= genome,
                          ...)
}

#' @describeIn vl_motif_counts Identify motifs in sequences
#' @export
vl_motif_counts.default <- function(sequences= NULL,
                                    pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds],
                                    genome,
                                    bg= "genome",
                                    p.cutoff= 5e-5)
{
  # Checks ----
  if(!"PWMatrixList" %in% class(pwm_log_odds))
    pwm_log_odds <- do.call(TFBSTools::PWMatrixList,
                            vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
  
  # Compute counts ----
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
  
  # Save ----
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
#' # Resize example peaks
#' SUHW <- vl_resizeBed(vl_SUHW_top_peaks, genome = "dm3")
#' STARR <- vl_resizeBed(vl_STARR_DSCP_top_peaks, genome = "dm3")
#' 
#' # Generate same number of random regions
#' random <- vl_control_regions_BSgenome(bed= STARR, genome= "dm3")
#' 
#' # Count JAPSPAR motifs (see below to use custom list of PWMs)
#' suhw <- vl_motif_counts(SUHW, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' starr <- vl_motif_counts(top_STARR, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' ctl <- vl_motif_counts(random, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' 
#' # Starting from sequence instead of bed file
#' seq <- vl_getSequence(SUHW, genome= "dm3")
#' seq_suhw <- vl_motif_counts(seq, genome= "dm3", pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' identical(suhw, seq_suhw)
#' 
#' # Motifs can also be counted using a custom PWMatrixList, for example for promoter motifs:
#' prom_db <- readRDS("/groups/stark/almeida/data/motifs/CP_motifs/CP_motifs_PWM.rds")
#' prom <- prom_db$Pwms_log_odds
#' for(i in seq(prom))
#'   prom[[i]]@profileMatrix <- vl_pwm_perc_to_log2(prom_db$Pwms_perc[[i]]@profileMatrix)
#' 
#' prom_motifs <- vl_motif_counts(STARR,
#'                                pwm_log_odds= prom,
#'                                genome= "dm3")
#' 
#' # Compute enrichment at SUHW peaks, using random controls as background
#' pl <- vl_motif_enrich(suhw,
#'                       ctl,
#'                       plot= F)
#' pl[order(padj)][1:3] # Top motif is su(Hw)
#' # Plot
#' plot(pl,
#'      top.enrich= 3)
#' 
#' 
#' # Compute enrichment at STARR peaks and plot at the same time
#' enr <- vl_motif_enrich(starr,
#'                        ctl,
#'                        plot= T,
#'                        order= "log2OR",
#'                        padj.cutoff= 1e-5)
#' # Positive enrichments identify typical S2 enhancer motifs
#' plot(enr[log2OR>.5])
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

#' Compute motif enrichment for several groups
#'
#' Compute motif enrichment for each element of a list (clusters, groups....)
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
#' # Resize example peaks
#' SUHW <- vl_resizeBed(vl_SUHW_top_peaks, genome = "dm3")
#' STARR <- vl_resizeBed(vl_STARR_DSCP_top_peaks, genome = "dm3")
#' 
#' # Generate same number of random regions
#' random <- vl_control_regions_BSgenome(bed= STARR, genome= "dm3")
#' 
#' # Combine the three sets of regions
#' combined <- rbindlist(list(SUHW= SUHW,
#'                            STARR= STARR,
#'                            random= random),
#'                       idcol = "cluster",
#'                       fill = T)
#' 
#' # Count JAPSPAR motifs
#' counts <- vl_motif_counts(combined, genome= "dm3", sel= pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds])
#' 
#' # Motifs can also be counted using a custom PWMatrixList, for example for promoter motifs:
#' prom_db <- readRDS("/groups/stark/almeida/data/motifs/CP_motifs/CP_motifs_PWM.rds")
#' prom <- prom_db$Pwms_log_odds
#' for(i in seq(prom)) # This set needs a fix
#'   prom[[i]]@profileMatrix <- vl_pwm_perc_to_log2(prom_db$Pwms_perc[[i]]@profileMatrix)
#' 
#' prom_motifs <- vl_motif_counts(combined,
#'                                pwm_log_odds= prom,
#'                                genome= "dm3")
#' 
#' # Enrichment per cluster
#' enr <- vl_motif_cl_enrich(split(counts, combined$cluster),
#'                           control.cl = "random")
#' 
#' # Plot
#' vl_par()
#' plot(enr,
#'      padj.cutoff= 1e-5)
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
#' @param pwms_log_odds A PWMatrixList (in log2 odd ratio format) for which motif matches should be mapped. Overrides sel and motifDB arguments (see above).
#' @param genome Genome to be used for coordinates ("dm6, "dm3") and as background for counting motifs when bg= "genome".
#' @param bg Background used to find motifs. Possible values include "genome" and "even". Default= "genome"
#' @param p.cutoff p.value cutoff used for motif detection. For enrichment analyses based on presence/absence of a motif, high cutoff might perform better (1e-4 or 5e-5) while for regression analyses, lower cutoffs might be prefered (5e-4). Default= 5e-5 (stringent).
#' @param collapse.overlapping Should overlapping motifs be merged? If TRUE (default), motif instances that overlap more than 70 percent of their width are collapsed.
#' 
#' @examples
#' # Find position of 3 different motifs within two regions
#' pos <- vl_motif_pos(vl_SUHW_top_peaks[1:2],
#'                     genome= "dm3",
#'                     pwm_log_odds= vl_Dmel_motifs_DB_full[c("cisbp__M2328", "flyfactorsurvey__suHw_FlyReg_FBgn0003567", "jaspar__MA0533.1"), pwms_log_odds, on= "motif_ID"])
#' 
#' # Starting from sequence
#' pos <- vl_motif_pos(sequence= "TGAGTTGTGTCTGAAATTGGGATTGCTGTTGCGACAATGCCTGTCTGACAGCATTGTCGATAAGAGCTTGAATCTGATTGGGGTCCATGGTAATATCTACCGTGGCACTATCTAACGGCCGACCTAATGCTTGGCCTACTTGCTCCTCCTCCCAGCTATCCTCGCTTTCGTATTCGACCTTAACCTTTCTGTAGTT#' ATGTGCCCAACTCATTGGTTGTTGGTTGGCACACCACAAATATACTGTTGCCGAGCACAATTGATCGGCTAAATGGTATGGCAAGAAAAGGTATGCAATATAATAATCTTTTATTGGGTATGCAACGAAAATTTGTTTCGTCAACGTATGCAATATTTTTTATTAAAAGAGGGTATGCAATGTATTTTATTAAAAACGGGTATGCAATATAATAATCTTTTATTGGG#' TATGCAACGAAAATTTGTTTCGTCAAAGTATGCAATATTTTTTATTAAAAGAGGGTATGCAATGTATTTTATTAAAAACGGGTATGCAATAAAAAATTATTTGGTTTCTCTAAAAAGTATGCAGCACTTATTTTTTGATAAGGTATGCAACAAAATTTTACTTTGCCGAAAATATGCAATGTTTTTGCGAATAAATTCAACGCACACTTATTACGTGGCCAGATACA#' CAACTTTTTTTTTTTTTTTTCACTCGTAAATTTCTTGATTGCGTCAAAGA",
#'                     genome= "dm3",
#'                     pwm_log_odds= vl_Dmel_motifs_DB_full[c("cisbp__M2328", "flyfactorsurvey__suHw_FlyReg_FBgn0003567", "jaspar__MA0533.1"), pwms_log_odds, on= "motif_ID"])
#' 
#' # Motifs can also be mapped using a custom PWMatrixList, for example for promoter motifs
#' prom_db <- readRDS("/groups/stark/almeida/data/motifs/CP_motifs/CP_motifs_PWM.rds")
#' prom <- prom_db$Pwms_log_odds
#' for(i in seq(prom))
#'   prom[[i]]@profileMatrix <- vl_pwm_perc_to_log2(prom_db$Pwms_perc[[i]]@profileMatrix)
#' 
#' # Select Ribo Protein gene promoters
#' proms <- rtracklayer::import("/groups/stark/annotations/dm3/dmel-all-filtered-r5.57_genes_and_transcripts_only.gff")
#' proms <- as.data.table(proms)
#' proms <- vl_resizeBed(proms[grepl("ribosomal protein", fullname)], "start", 150, 50)
#' proms[, seqnames:= paste0("chr", seqnames)]
#' 
#' prom_motifs <- vl_motif_pos(proms,
#'                             pwm_log_odds= prom,
#'                             genome= "dm3")
#' # Many instances of TCT motif are found, as expected
#' prom_motifs$TC_17_Zabidi
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
                                   pwm_log_odds= vl_Dmel_motifs_DB_full[collection=="jaspar", pwms_log_odds],
                                   genome,
                                   bg= "genome",
                                   p.cutoff= 5e-5,
                                   collapse.overlapping= TRUE)
{
  # Checks ----
  if(missing(genome))
    stop("genome is missing with no default")
  if(is.null(names(sequences)))
    names(sequences) <- seq(sequences)
  if(!"PWMatrixList" %in% class(pwm_log_odds))
    pwm_log_odds <- do.call(TFBSTools::PWMatrixList, pwms_log_odds)

  # Map motifs ----
  pos <- parallel::mclapply(pwm_log_odds,
                            function(x)
                            {
                              motifmatchr::matchMotifs(x,
                                                       sequences,
                                                       genome= genome,
                                                       p.cutoff= p.cutoff,
                                                       bg= bg,
                                                       out= "positions")[[1]]
                            },
                            mc.preschedule = FALSE,
                            mc.cores = data.table::getDTthreads()-1)

  # Name ----
  names(pos) <- TFBSTools::name(pwm_log_odds)
  
  # Format and collapse if specified ----
  pos <- lapply(pos, function(x)
  {
    x <- lapply(seq(x), function(i) {
      y <- as.data.table(x[[i]])
      y[, seqnames:= names(sequences)[i]]
      setcolorder(y, "seqnames")
      # Collapse motifs that overlap more than 70% (ignore strand)
      if(collapse.overlapping && nrow(y))
      {
        y <- vl_collapseBed(y,
                            min.gap = ceiling(mean(y$width)*0.7),
                            ignore.strand = TRUE)
        y <- y[, .(seqnames, start, end, width= end-start+1)]
      }
      return(y)
      })
    # Return
    names(x) <- names(sequences)
    return(x)
  })
  
  # Return ----
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
vl_pwm_perc_to_log2 <- function(perc_pwm)
{
  log2(perc_pwm/matrix(0.25, nrow= 4, ncol= ncol(perc_pwm))+1e-5)
}
