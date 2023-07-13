
#' Title
#'
#' @param file Path to bam file
#' @param sel columns to be imported (see Rsamtools::scanBamWhat()). Default= c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "mrnm", "mpos", "isize")
#' @param col.names Names of the columns to be imported. 
#'
#' @return Returns a bam file content as data.table
#' @export
#'
#' @examples
#' vl_importBam("path/to_bam_file")
vl_importBam <- function(file,
             sel= c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "mrnm", "mpos", "isize"),
             col.names= c("readID", "samFlag", "seqnames", "strand", "leftStart", "width", "mapq", "mateSeqnames", "mateLeftStart", "isize"))
{
  if(length(sel) != length(col.names))
  {
    col.names <- sel
    stop("Selected columns could not be renamed because length(sel) != length(col.names)")
  }
  param <- Rsamtools::ScanBamParam(what= sel)
  .c <- Rsamtools::scanBam(bam, param = param)[[1]]
  .c <- as.data.table(.c)
  setnames(.c, col.names)
  return(.c)
}