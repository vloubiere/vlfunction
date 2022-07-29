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

#' Dm6 GO database 
#'
#' Contains GO annotation from FB2020_05
#'
#' @usage Can be used with vl_GO_cluster. see ?vl_GO_cluster()
#' 
#' @format data.table containing GO metadata
#' @examples
#' GO: GO id
#' name: GO full name
#' type: GO type (biological_process, molecular_function, cellular_component)
#' FBgn: gene ID
#' Symbol: gene symbol
#' 
#' Object was generated using the following code:
#' tmp <- tempfile(fileext = ".gz")
#' download.file("http://ftp.flybase.net/releases/FB2020_05/precomputed_files/go/gene_association.fb.gz",
#' destfile = tmp)
#' go <- fread(tmp,
#' skip= 5,
#' select = c(2,3,5),
#' col.names = c("FBgn", "Symbol", "GO"))
#' tmp <- tempfile(fileext = ".gz")
#' download.file("http://ftp.flybase.net/releases/FB2020_05/precomputed_files/ontologies/go-basic.obo.gz",
#' destfile = tmp)
#' details <- ontologyIndex::get_ontology(tmp, extract_tags = "everything")
#' details <- data.table(GO= details$id,
#' name= details$name,
#' type= details$namespace)
#' vl_Dmel_GO_FB2020_05 <- merge(details, unique(go))
#' vl_Dmel_GO_FB2020_05 <- vl_Dmel_GO_FB2020_05[!(name %in% type)]
#' vl_Dmel_GO_FB2020_05[, type:= unlist(type)]
#' save(vl_Dmel_GO_FB2020_05,
#' file= "/groups/stark/vloubiere/vlfunction/data/vl_Dmel_GO_FB2020_05.RData")
"vl_Dmel_GO_FB2020_05"

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

#' Dm6 STRING database 
#'
#' Contains STRING annotations v11.5
#'
#' @usage Can be used with ?vl_STRING_interaction()
#' @examples
#' Code used to create DB
#' # Download STRING db v11.5 from Dmel
#' link <- "https://stringdb-static.org/download/protein.links.detailed.v11.5/7227.protein.links.detailed.v11.5.txt.gz"
#' tmp <- tempfile(fileext = ".txt.gz")
#' download.file(link, tmp)
#' DB <- fread(tmp)
#' # Download correspondance Dmel symbols
#' link <- "https://stringdb-static.org/download/protein.info.v11.5/7227.protein.info.v11.5.txt.gz"
#' download.file(link, tmp)
#' sym <- fread(tmp)
#' DB[sym, protein1_symbol:= i.preferred_name, on= "protein1==`#string_protein_id`"]
#' DB[sym, protein2_symbol:= i.preferred_name, on= "protein2==`#string_protein_id`"]
#' # Collapse unique interactions
#' DB[, merge:= paste0(sort(c(protein1_symbol, protein2_symbol)), collapse= ""), .(protein1_symbol, protein2_symbol)]
#' vl_Dmel_STRING_DB <- DB[, .SD[1], merge]
#' # SAVE
#' save(vl_Dmel_STRING_DB,
#' file = "/groups/stark/vloubiere/vlfunction/data/vl_Dmel_STRING_DB.RData", 
#' compress = "bzip2")
"vl_Dmel_STRING_DB"

#' mm10 STRING database 
#'
#' Contains STRING annotations v11.5
#'
#' @usage Can be used with ?vl_STRING_interaction()
#' @examples
#' Code used to create DB
#' # Download STRING db v11.5 from Dmel
#' link <- "https://stringdb-static.org/download/protein.links.detailed.v11.5/10090.protein.links.detailed.v11.5.txt.gz"
#' tmp <- tempfile(fileext = ".txt.gz")
#' download.file(link, tmp)
#' DB <- fread(tmp)
#' # Download correspondance Dmel symbols
#' link <- "https://stringdb-static.org/download/protein.info.v11.5/10090.protein.info.v11.5.txt.gz"
#' download.file(link, tmp)
#' sym <- fread(tmp)
#' DB[sym, protein1_symbol:= i.preferred_name, on= "protein1==`#string_protein_id`"]
#' DB[sym, protein2_symbol:= i.preferred_name, on= "protein2==`#string_protein_id`"]
#' # Collapse unique interactions
#' DB[, merge:= paste0(sort(c(protein1_symbol, protein2_symbol)), collapse= ""), .(protein1_symbol, protein2_symbol)]
#' vl_mm_STRING_DB <- DB[, .SD[1], merge]
#' # SAVE
#' save(vl_mm_STRING_DB,
#' file = "/groups/stark/vloubiere/vlfunction/data/vl_mm_STRING_DB.RData", 
#' compress = "bzip2")
"vl_mm_STRING_DB"
