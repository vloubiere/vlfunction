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
#' hit <- matchMotifs(vl_Dmel_motifs_DB$All_pwms_log_odds[selection], GRanges("chrX", IRanges(10000000, 10000500)), genome= "dm3", p.cutoff= 5e-4, bg= "even", out= "scores")
#' counts <- as.matrix(motifCounts(hit))
#' colnames(counts) <- name(vl_Dmel_motifs_DB$All_pwms_log_odds[selection])
#' counts <- as.data.table(counts)
#'
#' @source {"/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData"}
"vl_Dmel_motifs_DB"

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
#' @examples
#' hit <- matchMotifs(vl_CP_motifs_DB$Pwms_log_odds, GRanges("chrX", IRanges(10000000, 10000500)), genome= "dm3", p.cutoff= 5e-4, bg= "even", out= "scores")
#' counts <- as.matrix(motifCounts(hit))
#' colnames(counts) <- name(vl_CP_motifs_DB$Pwms_log_odds)
#' counts <- as.data.table(counts)
#'
#' @source {"/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData"}
"vl_CP_motifs_DB"


#' Dm6 GO database 
#'
#' Contains GO annotation from FB2020_05
#'
#' @usage Can be used with vl_GO_cluster. see ?vl_GO_cluster()
#' 
#' @format An object containing Motifs and related metadata
#' \describe{
#'   \item{GO}{GO ID}
#'   \item{name}{GO name}
#'   \item{type}{GO category}
#'   \item{partents}{parent GOs}
#'   \item{children}{children GOs}
#'   \item{Fbgn}{Gene flybase ID}
#'   \item{symbol}{Gene symbol}
#'   ...
#' }
#' 
#' @examples
#' Object was generated using the followin code:
#' go <- fread("../../genomes/dm6/FB2020_05_gene_association.fb",
#' skip= 5,
#' select = c(2,3,5))
#' colnames(go) <- c("FBgn", "Symbol", "GO")
#' go_details <- get_ontology("../../genomes/dm6/FB2020_05_go-basic.obo",
#' extract_tags = "everything")
#' go_details <- data.table(GO= go_details$id,
#' name= go_details$name,
#' type= go_details$namespace,
#' parents= go_details$parents,
#' children= go_details$children)
#' vl_fb_go_table_dm6_FB2020_05 <- go_details
#' save(vl_fb_go_table_dm6_FB2020_05, file= "../../vlfunction/data/vl_fb_go_table_dm6_FB2020_05.RData")
#' 
#' @source {"../../vlfunction/data/vl_fb_go_table_dm6_FB2020_05.RData"}
"vl_fb_go_table_dm6_FB2020_05"