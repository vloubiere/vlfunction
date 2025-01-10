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
  meta <- data.table(h5= h5,
                     bed= bed,
                     fa= fa)
  
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
  # Remove unused levels if any ----
  mot[, seqlvls:= factor(seqlvls, unique(seqlvls))]
  if(length(levels(mot$seqlvls))!=length(contrib.score))
    stop("The number of levels in seqlevels in mot should matche the length(contrib.score)")
  
  # Scale score from 0 to 1 and compute means ----
  min.score <- min(unlist(contrib.score), na.rm = TRUE)
  range.score <- diff(range(unlist(contrib.score), na.rm = TRUE))
  scaled.score <- lapply(contrib.score, function(x) (x-min.score)/range.score)
  mean.scaled.score <- sapply(scaled.score, mean)
  
  # Retrieve regions for rdm controls ----
  regions <- mot[, .(seqnames= levels(seqlvls), start= 1, end= seq.length)]
  mot.widths <- mot[!is.na(mot.count), round(mean(rbindlist(ir)$width)), motif]
  mot[mot.widths, bins.width:= i.V1, on= "motif"]
  
  # Loop over motifs and compute enrichment ----
  res <- mot[!is.na(bins.width), {
    # Bin all sequences
    bins <- vl_binBed(bed = regions,
                      bins.width = bins.width)
    # Select bins with correct width
    bins[, width:= end-start+1]
    bins <- bins[width==bins.width]
    bins[, seqlvls:= factor(seqnames, levels(seqlvls))]
    # Clean and add group
    bins$seqnames <- bins$binIDX <- NULL
    # For each motif
    .SD[, {
      # 2X random sampling
      sel <- sum(mot.count, na.rm= TRUE)*2
      set.seed(.GRP*sel)
      rdm <- bins[sample(.N, sel, replace = sel>.N)]
      # Motif instances
      mot <- .SD[!is.na(mot.count), ir[[1]], seqlvls]
      # Combine
      cmb <- list(motif= mot, rdm= rdm)
      cmb <- rbindlist(cmb, fill= TRUE, idcol = "group")
      # Extract motif scores
      cmb[, scores:= .(scaled.score[seqlvls])]
      cmb[, mot.scores:= .(.(scaled.score[[seqlvls]][start:end])), .(seqlvls, start, end)]
      # Compute enrichment
      cmb[, mean.score:= sapply(mot.scores, mean)]
      cmb[, log2OR:= log2(mean.score/mean.scaled.score[seqlvls])]
      # Compute OR and pval for each motif instance
      cmb[, pval:= {
        wilcox.test(mot.scores[[1]],
                    scores[[1]],
                    alternative= "greater")$p.value
      }, .(seqlvls, start, end)]
      cmb[, padj:= p.adjust(pval, method = "fdr")]
      # Fisher test
      .t <- table(motif= factor(cmb$group=="motif", c(TRUE, FALSE)),
                  signif= factor(cmb$padj<0.05, c(TRUE, FALSE)))
      OR.motif <- fisher.test(.t+1,
                              alternative = "greater")$estimate
      pval.motif <- fisher.test(.t,
                                alternative = "greater")$p.value
      # Select significant motifs coordinates
      sig.mot.coor <- cmb[group=="motif" & padj<0.05]
      sig.mot.coor <- sig.mot.coor[, .(seqlvls, start, end, width, mean.score, log2OR, padj)]
      # Return
      .(log2OR= log2(OR.motif),
        pval= pval.motif,
        sig.inst= .t[1,1],
        tot.inst= sum(.t[1,]),
        sig.rdm= .t[2,1],
        tot.rdm= sum(.t[2,]),
        sig.mot.coor= .(sig.mot.coor))
    }, motif]
  }, bins.width]
  # ADjusted p values
  res[, padj:= p.adjust(pval, "fdr")]
  
  # Clean ----
  clean <- res[, .(motif, log2OR, padj, 
                   sig.inst, tot.inst, sig.rdm, tot.rdm,
                   sig.mot.coor)]
  
  # Add missing values ----
  all <- data.table(motif= unique(mot$motif))
  final <- merge(all, clean, by= "motif", all.x= TRUE, sort= FALSE)
  
  # Return ----
  return(final)
}

#' Plot contribution scores matrix
#'
#' @param contrib.object A contribution object, as outputed by ?vl_importContrib()
#' @param enr And enrichment object, as outputed by ?vl_contrib_enrich() on the sequences of the contrib.object.
#' @param seqlvl A seqlvl present in in enr$sig.mot.coor$seqlvls.
#' @param min.count The minimum number of significant instances across (see enr$sig.inst). Default= 3L
#' @param best.by The group.by column for which only the instance with the best padjust should be returned.
#' @param sel If specified, only the motifs for which best.by %in% sel will be plotted.
#' @param xlab Default= "nt"
#' @param ylab Default= "Contribution"
#' @param xlim Default= sequence length.
#' @param ylim Default= range(contrib).
#'
#' @export
vl_plot_contrib_logo <- function(contrib.object,
                                 enr,
                                 seqlvl,
                                 min.count= 3L,
                                 best.by= "motif",
                                 sel,
                                 xlim,
                                 ylim,
                                 xlab= "nt",
                                 ylab= "Contribution")
{
  # Retrieve enriched motifs ----
  all.mots <- data.table::copy(enr)
  setnames(all.mots,
           old= best.by,
           new= "best.by")
  all.mots <- all.mots[, sig.mot.coor[[1]], .(best.by, sig.inst, log2OR, padj)]
  
  # Checks ----
  if(is.numeric(seqlvl) || length(seqlvl)>1)
    stop("seqlvl should be a character or a factor vector of length 1.")
  if(length(levels(all.mots$seqlvls)) != nrow(contrib.object))
    stop("The number of seqlvls in enr$sig.mot.coor should match the number of rows in contrib.object.")
  if(!seqlvl %in% levels(all.mots$seqlvls))
    stop("The seqlvl should be one of levels(enr$sig.mot.coor$seqlvls)")
  
  # Select sequence of interest ----
  mots <- all.mots[seqlvls==seqlvl & sig.inst>min.count, .SD[which.min(padj)], best.by]
  if(missing(sel))
    sel <- unique(mots$best.by)
  contrib <- contrib.object[which(levels(mots$seqlvls)==seqlvl)]
  contrib <- contrib[, .(base= unlist(tstrsplit(toupper(seq[[1]]), "")),
                         score= unlist(score))]
  
  # Plotting vars ----
  if(missing(xlim))
    xlim <- c(0, nrow(contrib))
  if(missing(ylim))
    ylim <- range(contrib$score)
  
  # Remove features outside xlim ----
  contrib[, xleft:= .I-1]
  contrib <- contrib[between(xleft, xlim[1], xlim[2])]
  mots <- mots[(start-1)<xlim[2] & end>xlim[1]]
  mots[(start-1)<xlim[1], start:= xlim[1]+1]
  mots[end>xlim[2], end:= xlim[2]]

  # Plotting ----
  plot(NA,
       xlim= xlim,
       ylim= ylim,
       xlab= xlab,
       ylab= ylab,
       frame= FALSE)
  contrib[, {
    vl_plotLetter(letter = base,
                  xleft = xleft,
                  ytop= score,
                  width = 1,
                  height = score)
  }, (contrib)]

  # Add motif boxes ----
  mots[best.by %in% sel, {
    if(.N)
    {
      rect(xleft = start-1,
           ybottom = ylim[1],
           xright = end,
           ytop = ylim[2])
      text((start-1+end)/2,
           ylim[2],
           best.by,
           pos= 3,
           xpd= T)
    }
  }]
  
  # Return motifs
  return(mots)
}


