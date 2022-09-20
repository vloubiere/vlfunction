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
#' par(mar= c(4,15,4,6))
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
#' @param add_motifs Show only n top enriched motifs
#' @param cex.width expansion factor for motif widths
#' @param cex.height expansion factor for motif heights
#' @return DT of enrichment values
#' @examples 
#' # Example run
#' suhw <- vl_motif_counts(vl_SUHW_top_peaks, "dm3")
#' ctls_regions <- vl_control_regions_BSgenome(bed= vl_SUHW_top_peaks, "dm3")
#' ctl <- vl_motif_counts(ctls_regions, "dm3")
#' par(mar= c(4,15,4,6))
#' pl <- vl_motif_enrich(suhw, ctl, add_motifs= T, padj_cutoff= 1e-3)
#' DT <- plot(pl, padj_cutoff= 1e-3)
#' vl_add_motifs(pl)
#' @return Returns a table of enrichment which can be plot using ?plot.vl_GO_enr
#' @export
vl_motif_enrich <- function(counts,
                            control_counts,
                            collapse_clusters= colnames(counts),
                            plot= T,
                            padj_cutoff= 0.05,
                            top_enrich= Inf, 
                            add_motifs= F,
                            cex.width= 1,
                            cex.height= 1)
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
  {
    DT <- plot(res,
               padj_cutoff= padj_cutoff,
               top_enrich= top_enrich)
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
#' @param add_motifs Show only n top enriched motifs
#' @param cex.width expansion factor for motif widths
#' @param cex.height expansion factor for motif heights
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
#' par(mar= c(4,20,4,4), las= 1)
#' pl <- vl_motif_cl_enrich(counts, 
#'                          cl_IDs = c(rep(1, nrow(top_SUHW)),
#'                                     rep(2, nrow(top_STARR))),
#'                          add_motifs= T)
#' pl <- plot(enr, padj_cutoff= 1e-5)
#' vl_add_motifs(pl$DT)
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
                               main= NA,
                               add_motifs= F,
                               cex.width= 1,
                               cex.height= 1)
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
    DT <- plot(res,
               padj_cutoff= padj_cutoff,
               top_enrich= top_enrich, 
               x_breaks= x_breaks,
               color_breaks= color_breaks,
               cex.balloons= cex.balloons,
               col= col,
               main= main)$DT
    if(add_motifs)
    {
      vl_add_motifs(DT,
                    cex.width= cex.width, 
                    cex.height= cex.height)
    }
  }
  
  return(res)
}

#' Add motifs to enrichment plots
#'
#' @param DT DT object output from vl_motif_enrich or vl_motif_cl_enrich
#' @param cex.width expansion factor for motif widths
#' @param cex.height expansion factor for motif heights
#' @export
vl_add_motifs <- function(DT, cex.width= 1, cex.height= 1)
{
  mats <- vl_Dmel_motifs_DB_full$pwms_perc[match(DT$motif_ID, vl_Dmel_motifs_DB_full$motif)]
  mats <- lapply(mats, TFBSTools::as.matrix)
  ax.width <- max(strwidth(DT$variable, cex= par("cex.axis")))+diff(grconvertX(c(0, par("mgp")[2]+0.5), "lines", "user"))
  coor <- vl_seqlogo(pwm = mats, 
                     x = par("usr")[1]-ax.width, 
                     y = DT$y, 
                     cex.width = cex.width,
                     cex.height = cex.height,
                     min_content= 0.05)
  coor[, segments(xleft, 
                  ybottom, 
                  xright, 
                  ybottom, 
                  xpd= T, 
                  lwd= 0.5)]
}

#' plot seqlogo
#'
#' plot seqlogo from pwm matrix
#'
#' @param pwm List of pwm matrices
#' @param x x positions
#' @param y positions (centered)
#' @param pos eith 2 (left) of 4 (right)
#' @param cex.width width expansion factor applied before plotting motifs
#' @param cex.height height expansion factor applied before plotting motifs
#' @param add Should the pwm be plot on the top of opened device? Default= T
#' @examples
#' pwm <- matrix(c(0.33,0.21,0.33,0.13,0.1,0.42,0.38,0.1,0.26,0.26,0.27,0.21,0,0.03,0.19,0.78,0.1,0.05,0.1,0.75,0.24,0.05,0.18,0.53,0.8,0.04,0,0.16,0.13,0.16,0.02,0.69,0.04,0.05,0.7,0.21,0.24,0.09,0.57,0.1,0.02,0.8,0.15,0.03,0.22,0.28,0.31,0.19,0.35,0.26,0.26,0.13,0.19,0.33,0.26,0.22), nrow= 4)
#' plot.new()
#' vl_seqlogo(pwm)
#' @export
vl_seqlogo <- function(pwm, 
                       x,
                       y,
                       pos= 2,
                       cex.width= 1,
                       cex.height= 1,
                       add= T,
                       min_content= 0)
{
  if(!pos %in% c(2,4))
    stop("Unsupported pos value. Use either 2 (left) or 4 (right)")
  if(is.matrix(pwm))
    pwm <- list(pwm)
  
  # Make object and index
  obj <- data.table(pwm, x, y, cex.width, cex.height)
  obj[, idx:= .I]
  
  # Width only depends on cex
  obj[, width:= strwidth("M", cex= cex.width), cex.width]
  
  # For each base, compute xleft, xright, ytop, ybottom
  pl <- obj[, {
    # Import PWM and melt
    .c <- as.data.table(pwm[[1]], keep.rownames= "base")
    .c <- melt(.c, id.vars = "base")
    # Compute motif content per column and normalize importance
    .c[, content:= sum(value*log2(value/c(0.25, 0.25, 0.25, 0.25)), na.rm= T), variable]
    .c[, norm:= value*(content/max(content))]
    # Remove flanks with little content (5% max)
    setorderv(.c, "variable")
    .c <- .c[min(which(norm>min_content)):max(which(norm>min_content))]
    # xleft depends on the pos (2 or 4)
    if(pos==4)
    {
      # Already correclty sorted earlier
      .c[, xleft:= x+((.GRP-1)*width), variable]
    }else if(pos==2)
    {
      setorderv(.c, "variable", -1)
      .c[, xleft:= x-(.GRP*width), variable]
    }
    # Rank from lowest to biggest importance -> inscreasing ytop pos
    setorderv(.c, "norm")
    .h <- strheight("M", cex= cex.height)
    .c[, c("height", "ytop"):= {
      heights <- norm*.h
      .(heights, (y-.h/2)+cumsum(heights))
    }, variable]
  }, .(idx, y, width)]

  # Plot
  pl[, vl_plotLetter(base[1], 
                     xleft[1], 
                     ytop[1], 
                     width[1], 
                     height[1]), .(base, xleft, ytop, width, height)]
  
  # Return object containing limits of each motif
  invisible(pl[, .(xleft= min(xleft), 
                   ybottom= min(ytop-height),
                   xright= max(xleft+width), 
                   ytop= max(ytop)), .(idx)])
}

#' plot seqlogo letter
#'
#' See function vl_seqlogo
#'
#' @param letter "A", "T", "C" or "G"
#' @param xleft xleft position
#' @param ytop ytop position
vl_plotLetter <- function(letter, xleft, ytop, width, height)
{
  if(letter=="T")
  {
    x <- c(0, 10, 10, 6, 6, 4, 4, 0) * 0.1
    y <- c(10, 10, 8.5, 8.5, 0, 0, 8.5, 8.5) * 0.1
    col <- "red"
  }else if(letter=="A")
  {
    x <- c(0, 4, 6, 2, 0, 4, 6, 10, 8, 4, 3.2, 6.8, 6.4, 3.6, 3.2) * 0.1
    y <- c(0, 10, 10, 0, 0, 10, 10, 0, 0, 10, 3, 3, 4, 4, 3) * 0.1
    col <- "forestgreen"
  }else if(letter=="C")
  {
    angle1 <- seq(0.3 + pi / 2, pi, length = 100)
    angle2 <- seq(pi, 1.5 * pi, length = 100)
    x.l1 <- 0.5 + 0.5 * sin(angle1)
    y.l1 <- 0.5 + 0.5 * cos(angle1)
    x.l2 <- 0.5 + 0.5 * sin(angle2)
    y.l2 <- 0.5 + 0.5 * cos(angle2)
    x.l <- c(x.l1, x.l2)
    y.l <- c(y.l1, y.l2)
    x <- c(x.l, rev(x.l))
    y <- c(y.l, 1 - rev(y.l))
    x.i1 <- 0.5 + 0.35 * sin(angle1)
    y.i1 <- 0.5 + 0.35 * cos(angle1)
    x.i1 <- x.i1[y.i1 <= max(y.l1)]
    y.i1 <- y.i1[y.i1 <= max(y.l1)]
    y.i1[1] <- max(y.l1)
    x.i2 <- 0.5 + 0.35 * sin(angle2)
    y.i2 <- 0.5 + 0.35 * cos(angle2)
    x.i <- c(x.i1, x.i2)
    y.i <- c(y.i1, y.i2)
    x1 <- c(x.i, rev(x.i))
    y1 <- c(y.i, 1 - rev(y.i))
    x <- c(x, rev(x1))
    y <- c(y, rev(y1))
    col <- "dodgerblue2"
  }else if(letter=="G")
  {
    angle1 <- seq(0.3 + pi / 2, pi, length = 100)
    angle2 <- seq(pi, 1.5 * pi, length = 100)
    x.l1 <- 0.5 + 0.5 * sin(angle1)
    y.l1 <- 0.5 + 0.5 * cos(angle1)
    x.l2 <- 0.5 + 0.5 * sin(angle2)
    y.l2 <- 0.5 + 0.5 * cos(angle2)
    x.l <- c(x.l1, x.l2)
    y.l <- c(y.l1, y.l2)
    x <- c(x.l, rev(x.l))
    y <- c(y.l, 1 - rev(y.l))
    x.i1 <- 0.5 + 0.35 * sin(angle1)
    y.i1 <- 0.5 + 0.35 * cos(angle1)
    x.i1 <- x.i1[y.i1 <= max(y.l1)]
    y.i1 <- y.i1[y.i1 <= max(y.l1)]
    y.i1[1] <- max(y.l1)
    x.i2 <- 0.5 + 0.35 * sin(angle2)
    y.i2 <- 0.5 + 0.35 * cos(angle2)
    x.i <- c(x.i1, x.i2)
    y.i <- c(y.i1, y.i2)
    x1 <- c(x.i, rev(x.i))
    y1 <- c(y.i, 1 - rev(y.i))
    x <- c(x, rev(x1))
    y <- c(y, rev(y1))
    h1 <- max(y.l1)
    r1 <- max(x.l1)
    h1 <- 0.4
    x.add <- c(r1, 0.5, 0.5, r1 - 0.2, r1 - 0.2, r1, r1)
    y.add <- c(h1, h1, h1 - 0.1, h1 - 0.1, 0, 0, h1)
    x <- c(rev(x), x.add)
    y <- c(rev(y), y.add)
    col <- "goldenrod1"
  }
  polygon(xleft+x*width, 
          ytop-(1-y)*height, 
          col= col,
          border= NA,
          xpd= T)
}

