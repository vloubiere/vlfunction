#' Merge bw files
#'
#' @param x fastq file path (read1)
#' @param output_folder Folder where to put processed files
#' @param output_prefix prefix for processed files
#' @examples 
#' vl_bw_merge(x= c("db/bw/cutnrun_reps/H3K27Ac_PH18_rep1.bw",
#' "db/bw/cutnrun_reps/H3K27Ac_PH18_rep2.bw"), 
#' output_folder = "db/bw/cutnrun_merge/",
#' output_prefix = "H3K27Ac_PH18_merge.bw")
#' @export
vl_bw_merge <- function(x, 
                        output_folder,
                        output_prefix)
{
  if(!is.character(x))
    stop("x should bw a character vector of bw paths")
  dat <- data.table(file= x)
  seqL <- list()
  dat <- dat[, {
    .c <- rtracklayer::import.bw(file)
    seqL[[.GRP]] <<- seqlengths(.c)
    as.data.table(.c)
  }, file]
  if(all(sapply(seqL, function(x) identical(x, seqL[[1]]))))
    seqL <- seqL[[1]] else
      stop("seqlenghts differ for some files!?!")
  res <- dat[, {
    .c <- matrix(sort(unique(c(start, end))), 
                 ncol = 2, 
                 byrow = T)
    colnames(.c) <- c("start", "end")
    as.data.table(.c)
  }, seqnames]
  res$score <- dat[res, sum(score), .EACHI, on= c("seqnames", "start<=end", "end>=start")]$V1
  res <- GRanges(res)
  seqlengths(res) <- seqL[match(seqlevels(res), names(seqL))]
  rtracklayer::export.bw(object = res, 
                         paste0(output_folder, "/", output_prefix, ".bw"))
}