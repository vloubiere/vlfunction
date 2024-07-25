#' Dmel motifs database
#'
#' Load vl_Dmel_motifs_DB_full (13,899 motifs) and vl_motifs_DB_v2 (6,502 motifs) DBs, assembled by Bernardo P. de almeida (https://github.com/bernardo-de-almeida/motif-clustering).
#' 
#' The vl_motifs_DB_v2 DB is a non-redundant version of vl_Dmel_motifs_DB_full.
#' 
#' @details
#' Columns:
#' \describe{
#'   \item{motif_name:}{Motif ID from TF_clusters_PWMs.RData}
#'   \item{FBgn:}{Curated FBgn symbols, separated by "/"}
#'   \item{motif_cluster:}{Motif cluster from Bernardo Almeida}
#'   \item{collection:}{Collection information (description not provided)}
#'   \item{collection_version:}{Collection version information (description not provided)}
#'   \item{species:}{Species name}
#'   \item{pwms_log_odds:}{List of PWM matrices expressed as log odds}
#'   \item{pwms_perc:}{List of PWM matrices expressed as percentages}
#' }
#' 
#' @usage
#' vl_Dmel_motifs_load()
#' Can be used with the \code{vl_motif_counts()} function to count motifs.
#'
#' @source {"/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData"}
#' @source 
#' \url{https://www.dropbox.com/scl/fi/w2p9ibid7xpcefpn5zvom/vl_Dmel_motifs_DB_full.RData?rlkey=764pj7stflznikbealoua2bvr&dl=1}
#' @source {"/groups/stark/almeida/Papers/DeepSTARR/Code/TF_motif_database/TF_clusters_PWMs.RData"}
#' @source
#' \url{https://www.dropbox.com/scl/fi/ut805p1vtb4l3y55gdo1k/vl_motifs_DB_v2.RData?rlkey=a4sm8rlv7t85zmmose8tbnocp&dl=1}
#' 
#' @export
vl_Dmel_motifs_load <- function()
{
  # Full dataset ----
  if(!exists("vl_Dmel_motifs_DB_full.RData"))
  {
    path <- paste0(tempdir(), "/vl_Dmel_motifs_DB_full.RData")
    if(!file.exists(path))
      download.file("https://www.dropbox.com/scl/fi/w2p9ibid7xpcefpn5zvom/vl_Dmel_motifs_DB_full.RData?rlkey=764pj7stflznikbealoua2bvr&dl=1",
                    path,
                    mode = "wb")
    load(path, envir = .GlobalEnv)
  }
  # Non redundant dataset ----
  if(!exists("vl_motifs_DB_v2.RData"))
  {
    path <- paste0(tempdir(), "/vl_motifs_DB_v2.RData")
    if(!file.exists(path))
      download.file("https://www.dropbox.com/scl/fi/ut805p1vtb4l3y55gdo1k/vl_motifs_DB_v2.RData?rlkey=a4sm8rlv7t85zmmose8tbnocp&dl=1",
                    path,
                    mode = "wb")
    load(path, envir = .GlobalEnv)
  }
}