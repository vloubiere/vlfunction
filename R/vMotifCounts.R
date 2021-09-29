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
                            resize= T,
                            extend= c(500, 500),
                            sel= "Dmel")
{
  # Checks
  if(is.data.table(bed))
    bed <- GRanges(bed)
  if(ncol(mcols(bed))>3)
    warning("provided bed file has many columns!! \n Given that the output can be massive, I would adive to reduce it to the minimu (coordinates + ID)")
  if(resize)
    if(length(extend)!=2 | !all(sign(extend)==1))
      stop("extend should be a vector containing two positive integers, e.g c(500, 500). If you want to use fancy resizing, resize before and set 'resize' argument to F!")
  
  # Select motifs
  if(is.null(sel))
    sel_IDs <- vl_Dmel_motifs_DB$metadata[,"motif_name"]
  if(identical("Dmel", sel))
  {
    idx <- !is.na(vl_Dmel_motifs_DB$metadata$Dmel) & # Associated to a known TF
      vl_Dmel_motifs_DB$metadata$X..motif_collection_name %in% # From a relevant DB
      c("flyfactorsurvey", "bergman", "jaspar", "idmmpmm", "cisbp")
    sel_IDs <- vl_Dmel_motifs_DB$metadata[idx, "motif_name"]
  }else if(all(sel %in% vl_Dmel_motifs_DB$metadata$motif_name))
    sel_IDs <- sel else
      stop("Sel should either be set to NULL (all motifs) or 'Dmel' (convenient set for Dmel) or should only contain IDs that exist in vl_Dmel_motifs_DB$metadata$motif_name")

  # Resize bed
  if(resize)
  {
    bed <- resize(bed, 1, "center")
    bed <- resize(bed, extend[1], fix = "end")
    bed <- resize(bed, extend[1]+extend[2], fix = "start")
  }
  
  # Compute counts
  sel_idx <- which(name(vl_Dmel_motifs_DB$All_pwms_log_odds) %in% sel_IDs)
  hit <- matchMotifs(vl_Dmel_motifs_DB$All_pwms_log_odds[sel_idx], 
                     bed, 
                     genome= genome, 
                     p.cutoff= 5e-4, 
                     bg= "even", 
                     out= "scores")
  counts <- as.matrix(motifCounts(hit))
  colnames(counts) <- name(vl_Dmel_motifs_DB$All_pwms_log_odds[sel_idx])
  counts <- as.data.table(counts)
  
  # Format output
  res <- cbind(as.data.table(bed), counts)
  res <- melt(res, 
              measure.vars = colnames(counts), 
              variable.name = "motif", 
              value.name = "motif_counts")
  res[as.data.table(vl_Dmel_motifs_DB$metadata), motif_name:= i.Dmel, on= "motif==motif_name"]
  res[, motif_name:= paste0(motif_name, "_", .SD[, rep(.GRP, .N), motif]$V1), motif_name]
  return(res)
}
  
