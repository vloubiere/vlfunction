#' Import contributions object
#'
#' @param h5 A vector of paths to h5 file(s) containing the contribution scores.
#' @param bed Bed file(s) containing the regions for which contributions were computed. By default, the file is searched in the same folder as the h5 file.
#' @param fa fasta file(s) containing the sequences for which contributions were computed. By default, the file is searched in the same folder as the h5 file.
#'
#' @return A contribution data.table containing, for each region, the contribution scores (as list) and the corresponding sequence.
#' @export
#'
#' @examples
#' dat <- vl_importContrib(h5= list.files("db/model_PH18/contributions/", "h5$", recursive = TRUE, full.names = TRUE))
#' 
vl_importContrib <- function(h5,
                             bed= list.files(dirname(h5), ".bed$", full.names = TRUE),
                             fa= list.files(dirname(h5), ".fa$", full.names = TRUE))
{
  # Metadata ----
  meta <- data.table(h5= h5, bed= bed, fa= fa)
  
  # Import ----
  dat <- meta[, {
    # bed regions
    .c <- vl_importBed(bed)
    if(!is.na(fa) && file.exists(fa))
      .c$seq <- as.character(seqinr::read.fasta(fa, as.string = TRUE)) else
        message("No valid fasta file found")
    if(!is.na(h5) && file.exists(h5))
    {
      .h <- rhdf5::h5read(h5, "contrib_scores/class")
      .c$score <- lapply(seq(dim(.h)[3]), function(i) rowSums(.h[,,i]))
    }else
      message("No valid h5 file found")
    .c
  }, .(h5, bed, fa)]
  dat$h5 <- dat$bed <- dat$fa <- NULL
  
  # Message ----
  if(any(nchar(dat$seq)>20))
  {
    options(datatable.prettyprint.char = 10)
    message("Printing option set to 10. To reset it to default, use options(datatable.prettyprint.char = NULL)")
  }

  # Return object ----
  return(dat)
}

#' Compute enrichment for each motif instance
#'
#' @param mot A list of motif created with !vl_motif_pos()
#' @param contrib.score A list of contribution scores, as in ?vl_importContrib() output
#' @param seq.length The length of the tiles used in the model. Default= 1001L.
#' @examples
#' # Contrib oject ----
#' files <- list.files(paste0("db/", model, "/contributions/"),
#'                     ".h5$",
#'                     recursive = T,
#'                     full.names = T)
#' contrib <- vl_importContrib(h5 = files)
#' 
#' # Find motifs ----
#' mot <- vl_motif_pos(sequences = contrib,
#'                     pwm_log_odds = vl_motifs_DB_v2[collection=="flyfactorsurvey" & !is.na(Dmel), pwms_log_odds],
#'                     bg = "genome",
#'                     genome = "dm6",
#'                     p.cutoff = 1e-4,
#'                     collapse.overlapping = TRUE)
#' 
#' # Motif enrichment ----
#' enr <- vl_contrib_enrich(mot = mot,
#'                          contrib.score = contrib$score,
#'                          seq.length = 1001L)
#' 
#' @return A data.table with the first columns corresponding to sequence levels (factor) and each column to the matches for the corresponding motif.
#' @export
vl_contrib_enrich <- function(mot,
                              contrib.score,
                              seq.length= 1001L)
{
  # Scale score from 0 to 1 and compute means ----
  min.score <- min(unlist(contrib.score), na.rm = TRUE)
  range.score <- diff(range(unlist(contrib.score), na.rm = TRUE))
  scaled.score <- lapply(contrib.score, function(x) (x-min.score)/range.score)
  mean.scaled.score <- sapply(scaled.score, mean)
  
  # Retrieve regions for rdm controls ----
  regions <- mot[, .(seqnames= seqlvls, start= 1, end= seq.length)]
  
  # Loop over motifs and compute enrichment ----
  motifs <- names(mot)[-1]
  res <- lapply(seq(motifs), function(i){
    # Import motif positions
    .c <- rbindlist(mot[, motifs[i], with= FALSE][[1]])
    .c[, width:= as.integer(width)]
    .c[, group:= "motif"]
    # Bin all sequences
    bins.width <- round(mean(.c$width))
    bins <- vl_binBed(bed = regions,
                      bins.width = bins.width)
    # Select bins with correct width
    bins[, width:= end-start+1]
    bins <- bins[width==bins.width]
    bins <- bins[end-start+1==bins.width]
    # 2X random sampling
    sel <- nrow(.c)*2
    set.seed(i)
    rdm <- bins[sample(.N, sel, replace = sel>.N)]
    rdm[, group:= "rdm"]
    # Combine
    cmb <- rbind(.c, rdm[, !"binIDX"])
    cmb[, seqlvls:= factor(seqnames, levels(mot$seqlvls))]
    # Extract motif scores
    cmb[, scores:= .(scaled.score[seqlvls])]
    cmb[, mot.scores:= .(.(scaled.score[[seqlvls]][start:end])), .(seqlvls, start, end)]
    # Compute enrichment
    cmb[, mean.mot.score:= sapply(mot.scores, mean)]
    cmb[, log2OR.inst:= log2(mean.mot.score/mean.scaled.score[seqlvls])]
    # Compute OR and pval for each motif instance
    cmb[, pval.inst:= {
      wilcox.test(mot.scores[[1]],
                  scores[[1]],
                  alternative= "greater")$p.value
    }, .(seqnames, start, end)]
    cmb[, padj.inst:= p.adjust(pval.inst, method = "fdr")]
    cmb$pval.inst <- NULL
    # Fisher test
    .t <- table(motif= factor(cmb$group=="motif", c(TRUE, FALSE)),
                signif= factor(cmb$padj.inst<0.05, c(TRUE, FALSE)))
    cmb$OR.motif <- fisher.test(.t+1,
                                alternative = "greater")$estimate
    cmb$pval.motif <- fisher.test(.t,
                                  alternative = "greater")$p.value
    cmb[, sig.inst:= .t[1,1]]
    cmb[, tot.inst:= sum(.t[1,])]
    cmb[, sig.rdm:= .t[2,1]]
    cmb[, tot.rdm:= sum(.t[2,])]
    
    # Return
    if(i==1 | i%%10==0)
      print(paste0(i, "/", length(motifs), " motifs processed."))
    return(cmb[group=="motif" & padj.inst<0.05])
  })
  
  # Compute padj per motif----
  names(res) <- motifs
  final <- rbindlist(res, idcol = "motif")
  final[, log2OR.motif:= log2(OR.motif)]
  padj <- unique(final[, .(motif, pval.motif)])
  padj[, padj.motif:= p.adjust(pval.motif, "fdr")]
  final[padj, padj.motif:= i.padj.motif, on= "motif"]
  
  # Clean and return ----
  clean <- final[, .(motif, log2OR.motif, padj.motif,
                     seqnames, start, end,
                     log2OR.inst, padj.inst,
                     sig.inst, tot.inst, sig.rdm, tot.rdm)]
  return(clean)
}

#' Plot contribution scores matrix
#'
#' @param bed A bed file containing a unique region for which contrib scores will be plotted.
#' @param h5 Path(s) to h5 files containing the contribution scores.
#' @param h5.bed Bed files containing the coordinates of the regions corresponding to provided h5 files.
#' @param genome The genome to be used.
#' @param agg.FUN In the case were several contribution scores would be found for a single nt, how should they be aggregated? Default= function(x) mean(x)
#' @param mot An optional bed file containing motifs to be added.
#' @param mot.name.column Name of the column containing the motif name.
#' @param xlab Default= "nt"
#' @param ylab Default= "Contribution"
#' @param xlim Default= sequence length
#' @param ylim Default= range(contrib)
#'
#' @export
#' 
#' @return contrib plot
#' 
vl_plot_contrib_logo <- function(bed,
                                 h5,
                                 h5.bed= list.files(dirname(h5), ".bed$", full.names = TRUE),
                                 h5.fa= list.files(dirname(h5), ".fa$", full.names = TRUE),
                                 genome,
                                 agg.FUN= function(x) mean(x),
                                 mot,
                                 mot.name.column= "motif_ID",
                                 xlab= "nt",
                                 ylab= "Contribution",
                                 xlim,
                                 ylim)
{
  # Import Bed
  bed <- vl_importBed(bed) # Very important!
  if(nrow(bed)>1)
    stop("Unique region should be specified")
  
  browser()
  # Metadata ----
  dat <- vl_importContrib(h5,
                          bed= h5.bed,
                          fa= h5.fa)
  
  # Intersect ----
  region <- dat[bed, on= c("seqnames", "start<=end", "end>=start")]
  
  
  # # Import contribution scores ----
  # dat <- meta[, {
  #   vl_importContrib(h5,
  #                    bed = h5.bed,
  #                    selection= bed)
  # }, .(h5, h5.bed)]
  # 
  # # Aggregate if necessary ----
  # if(uniqueN(dat[, .(seqnames, start)]) != nrow(dat))
  # {
  #   message("Some nucleotides had >1 contribution score assigned to it, which will be aggregated using agg.FUN")
  #   dat <- dat[, .(score= agg.FUN(score)), .(seqnames, start, end)]
  # }
  # 
  # # Get sequence ----
  # dat$base <- strsplit(vl_getSequence(bed, genome), "")[[1]]
  # 
  # # Plotting vars ----
  # dat[, xleft:= .I-1]
  # if(missing(xlim))
  #   xlim <- c(0, nrow(dat))
  # if(missing(ylim))
  #   ylim <- range(dat$score)
  # 
  # # Plotting ----
  # plot(NA,
  #      xlim= xlim,
  #      ylim= ylim,
  #      xlab= xlab,
  #      ylab= ylab,
  #      frame= FALSE)
  # 
  # dat[, {
  #   vl_plotLetter(base,
  #                 xleft = xleft,
  #                 ytop= score,
  #                 width = 1,
  #                 height = score)
  # }, (dat)]
  # 
  # # Add motif boxes ----
  # if(!missing(mot))
  # {
  #   mot <- vl_importBed(mot)
  #   mot <- vl_intersectBed(mot, bed, ignore.strand= TRUE)
  #   if(nrow(mot))
  #     mot[, {
  #       xl <- start-bed$start
  #       xr <- end-bed$start
  #       rect(xleft = xl,
  #            ybottom = ylim[1],
  #            xright = xr,
  #            ytop = ylim[2])
  #       text((xl+xr)/2,
  #            ylim[2],
  #            get(mot.name.column)[1],
  #            pos= 3,
  #            xpd= T)
  #     }, (mot)]
  # }
}


