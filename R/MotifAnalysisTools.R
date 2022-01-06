#' Counts motifs
#'
#' Counts motif occurences in a set of regions
#'
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

vl_motif_counts <- function(bed, 
                            genome= "dm6",
                            sel= vl_Dmel_motifs_DB_full[!is.na(vl_Dmel_motifs_DB_full$FBgn), motif])
{
  # Checks
  DT <- vl_importBed(bed)
  if(ncol(DT)>5)
    print("provided bed file has many columns!! \n Given that the output can be massive, I would advice to reduce it to the minimum (coordinates + ID)")
  if(any(!sel %in% vl_Dmel_motifs_DB_full$motif))
    stop("Some motif provided in 'sel' do not exist in vl_Dmel_motifs_DB_full$motif")

  # Select motifs
  sub <- vl_Dmel_motifs_DB_full[motif %in% sel]
  mot <- do.call(PWMatrixList, 
                 sub$pwms_log_odds)
  # Get sequence
  DT[, seq:= BSgenome::getSeq(BSgenome::getBSgenome(genome), seqnames, start, end, as.character= T)]
  res <- as.data.table(as.matrix(motifmatchr::matchMotifs(mot,
                                                          DT$seq,
                                                          p.cutoff= 5e-4,
                                                          bg= "even",
                                                          out= "scores")@assays@data[["motifCounts"]]))
  names(res) <- sub$motif
  res <- cbind(DT[, !"seq"], res)
  res <- melt(res,
              measure.vars = sub$motif,
              variable.name = "motif",
              value.name = "motif_counts")
  res <- merge(sub[, .(motif, motif_name, motif_FBgn= FBgn)],
               res)
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
#' @param auto_margin computes optimal plotting margins
#' 
#' @examples 
#' 
#' @return Fisher test data.table.
#' @export


vl_motif_cl_enrich_plot_only <- function(obj,
                                         padj_cutoff= 0.00001,
                                         log2OR_cutoff= 0,
                                         N_top= Inf,
                                         auto_margin= T)
{
  DT <- data.table::copy(obj)
  if(!all(c("motif", "motif_name", "cl", "pval", "padj", "log2OR") %in% names(DT)))
    stop("Input should contain c('motif', 'motif_name', 'cl', 'pval', 'padj', 'log2OR') columns. See ?vl_motif_cl_enrich() output")
  if(log2OR_cutoff<0)
    stop("log2OR_cutoff should be > 0")
  
  #----------------------------------#
  # Generate plot table
  #----------------------------------#
  # padj cutoff
  sel <- DT[, any(padj <= padj_cutoff & log2OR > log2OR_cutoff), motif][(V1), motif]
  if(length(sel)==0)
    stop("No enrichment found with provided paj cutoff!")
  # Plotting object
  pl <- DT[motif %in% sel & padj < 0.05 & log2OR > 0]
  pl[, '-log10(pval)':= -log10(pval)]
  # Select top motif/cluster
  setorderv(pl, "-log10(pval)", order = -1)
  top_motifs <- pl[, motif[seq(.N)<=N_top], cl]$V1
  pl <- pl[motif %in% top_motifs]
  # Handle infinite log2OR
  if(any(is.infinite(pl$log2OR)))
  {
    warning("There are suspcious infinite log2OR values found after padj cutoff")
    print(unique(pl[is.infinite(pl$log2OR), motif]))
  }
  pl[log2OR==Inf, log2OR:= max(pl[is.finite(log2OR), log2OR])]
  pl[log2OR==(-Inf), log2OR:= min(pl[is.finite(log2OR), log2OR])]
  # pval coloring
  pval_lims <- range(pl$`-log10(pval)`)
  Cc <- colorRamp2(pval_lims, 
                   colors = c("blue", "red"))
  pl[, col:= Cc(`-log10(pval)`)]
  # Points scaling factor
  pl[, size:= (2^log2OR)]
  pl[, size:= par("cex")*log2(size+1)]
  size_lims <- range(pl$size)
  # Y coordinates
  setorderv(pl, 
            c("cl", "-log10(pval)", "log2OR", "motif"), 
            order = c(1, -1, -1, 1))
  pl[, y:= as.numeric(min(.I)), motif]
  pl[, y:= seq(1, 0, length.out = length(unique(pl$motif)))[.GRP], keyby= y]
  # X coordinates
  pl[, x:= seq(0, 1, length.out = length(unique(pl$cl)))[.GRP], keyby= cl]
  
  #-----------------------------#
  # PLOT
  #-----------------------------#
  if(auto_margin)
    par(mai = c(max(strwidth(pl$cl, "inches"))+0.5,
                max(strwidth(pl$motif_name, "inches"))+0.5,
                0.5,
                strwidth("-log10(pval)", "inches")+0.25),
        xaxs= "i",
        yaxs= "i")
  
  # Lines
  plot.new()
  segments(unique(pl$x),
           0,
           unique(pl$x),
           1,
           xpd=T)
  segments(0,
           unique(pl$y),
           1,
           unique(pl$y),
           xpd=T)
  # Points
  points(pl$x, 
         pl$y, 
         cex= pl$size, 
         col= adjustcolor(pl$col, 0.8),
         pch= 16,
         xpd= T)
  axis(1,
       lty= 0,
       at = unique(pl$x), 
       labels = unique(pl$cl), 
       las= 2)
  axis(2,
       lty= 0,
       at = unique(pl[, .(motif_name, y)])$y, 
       labels = unique(pl[, .(motif_name, y)])$motif_name,
       las= 2)
  # Legend pval
  xleft <- grconvertX(1, "npc", "inches")+grconvertX(1, "lines", "inches")
  xright <- xleft+grconvertX(1, "lines", "inches")
  xleft <- grconvertX(xleft, "inches", "npc")
  xright <- grconvertX(xright, "inches", "npc")
  ybottom <- 0.7
  ytop <- grconvertY(1, "npc", "inches")-grconvertY(1, "chars", "inches")
  ytop <- grconvertY(ytop, "inches", "npc")
  rasterImage(matrix(rev(Cc(seq(min(pval_lims), max(pval_lims), length.out = 101)))),
              xleft,
              ybottom,
              xright,
              ytop,
              xpd=T)
  ticks <- axisTicks(pval_lims, log=F)
  ymin.ticks <- ybottom+(min(ticks)-pval_lims[1])/diff(pval_lims)*(ytop-ybottom)
  ymax.ticks <- ybottom+(max(ticks)-pval_lims[1])/diff(pval_lims)*(ytop-ybottom)
  text(xright,
       seq(ymin.ticks, ymax.ticks, length.out = length(ticks)),
       labels = ticks,
       pos=4,
       xpd= T,
       cex= 0.6,
       offset= 0.25)
  text(xleft,
       1-0.5*strheight("A"),
       labels = "-log10(pval)",
       pos=4,
       xpd= T,
       cex= 0.8,
       offset= 0)
  
  # Legend balloons
  scale <- axisTicks(size_lims, log= F)
  maxBalloonInch <- strheight("A", units = "inches", cex= max(scale))*0.75
  bx <- grconvertX(xleft, "npc", "inches")+maxBalloonInch/2
  bx <- grconvertX(bx, "inches", "npc")
  btop <- grconvertY(0.6, "npc", "inches")-maxBalloonInch/2
  bbot <- btop-maxBalloonInch*(length(scale)-1)
  btop <- grconvertY(btop, "inches", "npc")
  bbot <- grconvertY(bbot, "inches", "npc")
  points(rep(bx, length(scale)),
         seq(bbot, btop, length.out = length(scale)), 
         xpd= T,
         col= "black",
         cex= scale,
         pch= 16)
  bx <- grconvertX(bx, "npc", "inches")+maxBalloonInch/2
  bx <- grconvertX(bx, "inches", "npc")
  text(rep(bx, length(scale)),
       seq(bbot, btop, length.out = length(scale)),
       labels= scale,
       pos= 4,
       xpd= T,
       offset= 0,
       cex= 0.8)
  tity <- grconvertY(0.6, "npc", "inches")+grconvertY(1, "chars", "inches")/2
  tity <- grconvertY(tity, "inches", "npc")
  text(xleft,
       tity,
       labels = "log2OR",
       pos= 4,
       xpd= T,
       cex= 0.8,
       offset= 0)
  
  invisible(pl)
}

#' Compute motif enrichment
#'
#' Compute motif enrichment for the cluster in cl_columns, using all the lines as background
#'
#' @param obj A motif count object similar to the output of ?vl_motif_counts() (colnames:'motif', 'motif_counts', 'motif_name') + one cluster column! 
#' @param cl_column name of the cluster column
#' @param bg cluster IDs to be used as background. Default unique(obj[[cl_column]])
#' @param comp_expr Expression used for contingency table. Default "motif_count>0"
#' 
#' @examples 
#' 
#' @return Fisher test data.table.
#' @export

vl_motif_cl_enrich <- function(obj, 
                               cl_column,
                               bg= unique(obj[[cl_column]]),
                               comp_expr= "motif_counts>0")
{
  DT <- data.table::copy(obj)
  # Checks
  if(!cl_column %in% names(DT))
    stop("cl_column does not exist in provided obj") else if(cl_column != "cl")
      names(DT)[names(DT)==cl_column] <- "cl"
  if(class(DT$cl) != class(bg))
    stop("cluster column and bg are not from same class -- > coerce?")
  if(!all(c("motif", "motif_name", "motif_FBgn", "motif_counts") %in% names(DT)))
    stop("Provided DT should contain c('motif', 'motif_name', 'motif_FBgn', 'motif_counts') columns. see ?vl_motif_counts() output!")
  if(anyNA(bg))
    stop("bg contains NAs. cl_column should not contain NAs OR they should not be used for comparisons!")
  DT <- DT[, .(motif, motif_name, motif_counts, motif_FBgn, cl)]

  #-----------------------#
  # Format result table
  #-----------------------#
  cmb <- DT[, {
    data.table(V1= na.omit(unique(cl))) # make DT containing clusters to test
  }, .(mot= motif, motif_name, motif_FBgn)]
  
  #-----------------------#
  # For each motif/cluster combination, compute association using fisher
  #-----------------------#
  cmb[, c("OR", "pval"):= {
    # Extract motif from DT, restrict to regions from tested cluster OR bg, cast contingency table
    .t <- table(DT[motif==mot & cl %in% c(V1, bg), .(cl==V1, 
                                                     eval(parse(text= comp_expr)))])
    # If motif present in both tested cl and bg, do fisher test
    if(identical(dim(.t), as.integer(c(2,2))))
      res <- as.data.table(as.list(fisher.test(.t)[c("estimate", "p.value")])) else
        data.table(numeric(), numeric())
  }, .(mot, V1)]
  cmb <- na.omit(cmb)
  cmb[, padj:= p.adjust(pval, method = "fdr"), pval]
  cmb[, log2OR:= log2(OR)]
  setnames(cmb, 
           c("mot", "V1"),
           c("motif", "cl"))
  return(cmb)
}
