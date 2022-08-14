#' Dmel motifs data base 
#'
#' Contains 13899 motifs usable for motif enrichments in Dmel
#'
#' @usage Can be used with the ?vl_motif_counts() function to count motifs
#' 
#' @format An object containing 13899 motifs and related metadata
#' \describe{
#'   \item{motif}{motif IDs consitent with Bernardo's ones.}
#'   \item{motif_name}{Curated Dmel symbols. sep= "/"}
#'   \item{FBgn}{Curated FBgn symbols. sep= "/"}
#'   \item{collection}{Motif collection.}
#'   \item{Motif_cluster_name}{Cluster names with the respective motif types.}
#'   \item{Pwms_log_odds}{PWM expressed as log odds}
#'   \item{Pwms_perc}{PWM expressed as percentage}
#' }
#' @source {"/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData"}
"vl_Dmel_motifs_DB_full"

#' SUHW top peaks
#'
#' Example set that can be used for motif enrichment...
#'
#' @usage see ?vl_motif_enrich()
#' @format Narrowpeaks file containing top SUHW peaks on dm3 chanonical chromosomes
"vl_SUHW_top_peaks"

#' STARR-Seq top peaks
#'
#' Example set that can be used for motif enrichment...
#'
#' @usage see ?vl_motif_enrich()
#' @format Narrowpeaks file containing top STARR-Seq peaks (DSCP) on dm3
"vl_STARR_DSCP_top_peaks"

#' Set of genes
#'
#' Example set that can be used for GO enrichment. Contain RpL and HOX genes
#'
#' @usage see ?vl_GO_enrich()
#' @format FBgn, symbol and GO (RpL/HOX)
"vl_genes_set"

#' Thermofisher enzymes
#'
#' Contains cutting sites for most commercial enzymes
#'
#' @usage Can be used with vl_digest. see ?vl_digest()
#' 
#' @format A table containing thermofisher enzymes and related info
#' \describe{
#'   \item{name1}{enzyme name}
#'   \item{name2}{Alternative names}
#'   \item{fastdigest}{FastDigest?}
#'   \item{buffer}{Buffer to be used}
#'   \item{temperature}{Temeprature for digest}
#'   \item{cutsite}{cutsite pattern}
#'   \item{consensus_F}{consensu motif}
#'   ...
#' }
#' @examples 
#' The table was generated using the following code:
#' setwd("/groups/stark/vloubiere/projects/z_miscellaneous/thermofisher_restriction_enzymes_table/")
#' require(rvest)
#' html <- read_html("https://www.thermofisher.com/at/en/home/brands/thermo-scientific/molecular-biology/thermo-scientific-restriction-modifying-enzymes/restriction-enzymes-thermo-scientific/conventional-restriction#' #' #' -enzymes-thermo-scientific.html")
#' dat <- as.data.table(matrix(html %>%
#'                               html_nodes("#quickchart a") %>%
#'                               html_text(), ncol= 3, byrow = T))
#' colnames(dat) <- c("name1", "name2", "fastdigest")
#' dat[, buffer := html %>%
#'       html_nodes(".optimalbuffer") %>%
#'       html_text()]
#' dat[, temperature:= html %>%
#'       html_nodes(".otherconditions tr:nth-child(1) td") %>%
#'       html_text()]
#' dat[, cutsite := html %>%
#' html_nodes("code") %>%
#'       html_text()]
#' dat[, consensus_F := {
#'   .c <- unlist(strsplit(cutsite, ""))
#'   .c <- .c[.c %in% DNA_ALPHABET[1:15]]
#'   paste0(.c, collapse= "")}, cutsite]
#' @source {"https://www.thermofisher.com/at/en/home/brands/thermo-scientific/molecular-biology/thermo-scientific-restriction-modifying-enzymes/restriction-enzymes-thermo-scientific/conventional-restriction-enzymes-thermo-scientific.html"}
"vl_thermofisher_restriction_enzymes_table"