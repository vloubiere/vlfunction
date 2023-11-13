#' Plot enhancer with its motifs
#'
#' @param sequences Named character vector of sequences to analyse. If provided takes over bed argument (in the case where both are specified)
#' @param bed Either a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param genome Genome to be used for coordinates ("dm6, "dm3")
#' @param p.cutoff Pval cutoff used for motif detection
#' @param sel List of motifs to compute. see vl_Dmel_motifs_DB_full$motif
#' @param xleft Left position for plotting
#' @param ybottom Bottom position for plotting
#' @param width Plot width
#' @param height Plot height
#' @param col Color (1 per motif)
#' 
#' @examples
#' plot.new()
#' vl_plot_enh_motifs(bed= vl_SUHW_top_peaks[1:2], 
#'                    genome= "dm3", 
#'                    sel= c("cisbp__M2328", "flyfactorsurvey__suHw_FlyReg_FBgn0003567", "jaspar__MA0533.1"), 
#'                    xleft=0, 
#'                    ybottom=c(0.25,0.5), 
#'                    width=1, 
#'                    height= 0.1)
#' 
#' @export
vl_plot_enh_motifs <- function(sequences, ...) UseMethod("vl_plot_enh_motifs")

#' @describeIn vl_plot_enh_motifs Plot enhancer with its motifs
#' @export
vl_plot_enh_motifs.data.table <- function(bed, genome, ...)
{
  sequences <- vl_getSequence(bed, genome)
  vl_plot_enh_motifs.character(sequences, ...)
}

#' @describeIn vl_plot_enh_motifs Plot enhancer with its motifs
#' @export
vl_plot_enh_motifs.character <- function(sequences,
                                         p.cutoff= 5e-4,
                                         sel,
                                         xleft,
                                         ybottom,
                                         width,
                                         height,
                                         col= rainbow(length(sel)),
                                         col_alpha= 0.3)
{
  if(any(!sel %in% vl_Dmel_motifs_DB_full$motif_ID))
    stop("Some motif_ID provided in 'sel' do not exist in vl_Dmel_motifs_DB_full$motif_ID")
  sub <- vl_Dmel_motifs_DB_full[match(unique(sel), motif_ID)]
  sub[, col:= rep(col, length.out= .NGRP)[.GRP], motif_ID]
  sub <- sub[, .(sequence= seq(length(sequences))), .(motif_ID, col)]
  sub[, length:= nchar(sequences)[.GRP], sequence]
  sub[, xleft:= rep(xleft, length.out= .NGRP)[.GRP], sequence]
  sub[, ybottom:= rep(ybottom, length.out= .NGRP)[.GRP], sequence]
  sub[, width:= rep(width, length.out= .NGRP)[.GRP], sequence]
  sub[, height:= rep(height, length.out= .NGRP)[.GRP], sequence]
  
  # Get positions
  pos <- vl_motif_pos(sequences, p.cutoff= p.cutoff, sel= sel)

  # Plot
  pl <- sub[, unlist(pos, recursive = F)[[.GRP]][, !"width"], (sub)]
  pl[, rect(xleft+start/length*width,
            ybottom,
            xleft+end/length*width,
            ybottom+height,
            col= adjustcolor(col, col_alpha),
            border= NA)]
  rect(xleft,
       ybottom,
       xleft+width,
       ybottom+height,
       border= "grey20")

  invisible(pl)
}

