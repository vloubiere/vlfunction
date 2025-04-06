#' Import samtools using rsamtools
#'
#' @param file Path to bam file
#' @param sel columns to be imported (see Rsamtools::scanBamWhat()). Default= c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "mrnm", "mpos", "isize")
#' @param col.names Names of the columns to be imported. Default= c("readID", "samFlag", "seqnames", "strand", "leftStart", "width", "mapq", "mateSeqnames", "mateLeftStart", "isize")
#'
#' @return Returns a bam file content as data.table
#' @export
#'
#' @examples
#' importBam("path/to_bam_file")
importBam <- function(file,
                      sel= c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "mrnm", "mpos", "isize"),
                      col.names= c("readID", "samFlag", "seqnames", "strand", "leftStart", "width", "mapq", "mateSeqnames", "mateLeftStart", "isize"))
{
  # Checks
  if(length(sel) != length(col.names)) {
    col.names <- sel
    stop("Selected columns could not be renamed because length(sel) != length(col.names)")
  }

  # Import file
  param <- Rsamtools::ScanBamParam(what= sel)
  .c <- Rsamtools::scanBam(file, param = param)[[1]]

  # Coerce special elements (typically DNAStringSet / PhredQuality) to character
  toCoerce <- which(!sapply(.c, function(x) is.vector(x) | is.factor(x)))
  for(i in toCoerce){
    .c[[i]] <- as.character(.c[[i]])
  }

  # To data.table
  .c <- as.data.table(.c)
  setcolorder(.c, sel)
  setnames(.c, old = sel, new= col.names)

  # Return
  return(.c)
}
