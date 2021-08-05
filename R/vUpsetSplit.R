#' upset Split
#'
#' Split a list of objects in upset-plot lie classes
#'
#' @param dat_list List containing the names to intersect (see examples)
#' @examples 
#' test <- list(A= 1:1000, B= 1:1500, C= 1000:1750)
#' vl_upset_split(test)
#' @export

vl_upset_split <- function(dat_list)
{
  if(is.null(names(dat_list)))
    stop("The list should be named!")
  if(any(grepl("\\|", names(dat_list))))
    stop("list names should not contain any '|' cause they are use internally")
  if(max(nchar(names(dat_list)))<3)
    names(dat_list) <- paste0("  ", names(dat_list), "  ")
  # if(par("mfrow"))

  #-------------------------#
  # Main barplot
  #-------------------------#
  dat <- rbindlist(lapply(dat_list, as.data.table), idcol = T)
  colnames(dat)[2] <- "intersect"
  dat <- dat[, .(.id= paste(.id, collapse = "|")), intersect]
  dat <- dat[, .(N= .N, elements= .(intersect)), .id]
  return(dat)
}
