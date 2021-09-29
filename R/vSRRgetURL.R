#' Extract SRR download fastq
#'
#' Tries to extract the download link of the SRA file for given SRR
#'
#' @param SRR SRR ID
#' @examples 
#' vl_SRR_url("SRR13325541") 
#' @return A character vector containing the SRA download link

#' @export

vl_SRR_url <- function(SRR)
{
  if(!"rvest" %in% rownames(installed.packages()))
    stop("Install rvest package")
  .c <- rvest::read_html(paste0("https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=", SRR))
  url <- .c %>%
    rvest::html_nodes("#sra-viewer-app .first a") %>%
    rvest::html_text()
  return(url)
}
