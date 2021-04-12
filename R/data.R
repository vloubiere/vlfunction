#' Dmel motifs data base 
#'
#' Contains 13899 motifs usable for motif enrichments in Dmel
#'
#' @usage Can be used with the matchMotifs function from the motifmatchr package (see example)
#' 
#' @format An object containing 13899 motifs and related metadata
#' \describe{
#'   \item{Motif_cluster}{901 clusters containing redunsant motifs. Used collections: bergman, cisbp, flyfactorsurvey, homer, jaspar, stark, idmmpmm.}
#'   \item{Motif_cluster_name}{Cluster names with the respective motif types.}
#'   \item{Pwms_log_odds}{PWM expressed as log odds}
#'   \item{Pwms_perc}{PWM expressed as percentage}
#' }
#' 
#' @examples
#' selection <- 1:100
#' hit <- matchMotifs(Dmel_motifs_DB$All_pwms_log_odds[selection], GRanges("chrX", IRanges(10000000, 10000500)), genome= "dm3", p.cutoff= 5e-4, bg= "even", out= "scores")
#' counts <- as.matrix(motifCounts(hit))
#' colnames(counts) <- name(Dmel_motifs_DB$All_pwms_log_odds[selection])
#' counts <- as.data.table(counts)
#'
#' @source {"/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData"}
"Dmel_motifs_DB"

#' Core Promoter motifs data base 
#'
#' Contains 19 core promoter motifs (Vanja + Bernie)
#'
#' @usage Can be used with the matchMotifs function from the motifmatchr package (see example)
#' 
#' @format An object containing Motifs and related metadata
#' \describe{
#'   \item{Pwms_log_odds}{PWM expressed as log odds}
#'   \item{Pwms_perc}{PWM expressed as percentage}
#'   ...
#' }
#' 
#' #' @examples
#' hit <- matchMotifs(CP_motifs_DB$Pwms_log_odds, GRanges("chrX", IRanges(10000000, 10000500)), genome= "dm3", p.cutoff= 5e-4, bg= "even", out= "scores")
#' counts <- as.matrix(motifCounts(hit))
#' colnames(counts) <- name(CP_motifs_DB$Pwms_log_odds)
#' counts <- as.data.table(counts)
#'
#' @source {"/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData"}
"CP_motifs_DB"