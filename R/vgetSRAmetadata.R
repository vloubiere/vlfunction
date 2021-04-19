#' get SRA metadata 
#'
#' Donwload SRA metadata info from GSE
#'
#' @param GSE GSE number of interest. 
#' @examples
#' vl_sra_metadata(GSE= "GSE119708", SRAdb= "/groups/stark/vloubiere/exp_data/SRAmetadb.sqlite")
#' @return metadata data.table object
#' @format metadata data.table object
#' \describe{
#'   \item{study_name}{GSE id}
#'   \item{run}{SRR ID}
#'   \item{library_layout}{"PAIRED - " or "SINGLE - "}
#'   \item{experiment_title}{experiment_title}
#'   \item{ftp}{List of ftp download link}
#' }
#' @export

vl_sra_metadata <- function(GSE,
                            SRAdb= "/groups/stark/vloubiere/exp_data/SRAmetadb.sqlite")
{
  if(!file.exists(SRAdb))
    stop("SRAdb file does not exist. See ?getSRAdbFile")
  con <- dbConnect(SQLite(), SRAdb)
  .c <- as.data.table(getSRA(GSE, sra_con = con))
  .c[, ftp:= .(list(as.character(getFASTQinfo(in_acc = run, srcType = "ftp", sra_con = con)$ftp))), (.c)]
  return(.c)
}
