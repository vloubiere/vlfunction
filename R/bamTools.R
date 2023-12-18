
#' Title
#'
#' @param file Path to bam file
#' @param sel columns to be imported (see Rsamtools::scanBamWhat()). Default= c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "mrnm", "mpos", "isize")
#' @param col.names Names of the columns to be imported. Default= c("readID", "samFlag", "seqnames", "strand", "leftStart", "width", "mapq", "mateSeqnames", "mateLeftStart", "isize")
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
  .c <- Rsamtools::scanBam(file, param = param)[[1]]
  if(any(sapply(.c, class)=="DNAStringSet"))
  {
    idx <- which(sapply(.c, class)=="DNAStringSet")
    for(i in idx)
      .c[[i]] <- as.character(.c[[i]])
  }
  .c <- as.data.table(.c)
  setnames(.c, col.names)
  return(.c)
}

#' Title
#'
#' Uses samtools to import a bam file
#' @param file bam file path
#' @param extra_arg Extra arg to be passed to samtools view
#'
#' @return Imported bam file
#' @export
#'
#' @examples
#' vl_importBamRaw("path/to/bam/file.bam")
vl_importBamRaw <- function(file,
                            extra_arg)
{
  cmd <- "module load build-env/2020; module load samtools/1.9-foss-2018b; samtools view"
  if(!missing(extra_arg))
    cmd <- paste(cmd,
                 extra_arg)
  cmd <- paste(cmd, file)
  fread(cmd= cmd,
        fill= T)
}
