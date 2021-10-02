#' Counts motifs
#'
#' Counts motif occurences in a set of regions
#'
#' @param bed Either a GRanges object or a data.table that can be coerced to it
#' @param genome Genome to be used for coordinates ("dm6, "dm3")
#' @param resize Should the regions be resize according to extend arg? Default= T
#' @param extend Vector containing two positive integers (e.g c(500,500) indicating how the regions should be extend around their center
#' @param sel Either "Dmel" (convenient set of Dmel motifs), NULL (all motifs) or a list of IDs existing in vl_Dmel_motifs_DB$metadata$motif_name
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
  if(!is.data.table(bed))
    bed <- as.data.table(bed)
  if(ncol(bed)>5)
    print("provided bed file has many columns!! \n Given that the output can be massive, I would advice to reduce it to the minimum (coordinates + ID)")
  if(any(!sel %in% vl_Dmel_motifs_DB_full$motif))
    stop("Some motif provided in 'sel' do not exist in vl_Dmel_motifs_DB_full$motif")
  
  # Select motifs
  sub <- vl_Dmel_motifs_DB_full[motif %in% sel]
  mot <- do.call(PWMatrixList, 
                 sub$pwms_log_odds)
  # Get sequence
  bed[, seq:= getSeq(getBSgenome(genome), seqnames, start, end, as.character= T)]
  res <- as.data.table(as.matrix(matchMotifs(mot,
                                             bed$seq,
                                             p.cutoff= 5e-4,
                                             bg= "even",
                                             out= "scores")@assays@data[["motifCounts"]]))
  names(res) <- sub$motif
  res <- cbind(bed[, !"seq"], res)
  res <- melt(res,
              measure.vars = sub$motif,
              variable.name = "motif",
              value.name = "motif_counts")
  res <- merge(sub[, .(motif, motif_name)],
               res)
  return(res)
}

