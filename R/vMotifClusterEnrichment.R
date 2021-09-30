#' Compute motif enrichment
#'
#' Compute motif enrichment for the cluster in cl_columns, using all the lines as background
#'
#' @param obj A motif count object similar to the output of ?vl_motif_counts() (colnames:'motif', 'motif_counts', 'motif_name') + one cluster column! 
#' @param cl_column name of the cluster column
#' 
#' @examples 
#' 
#' @return Fisher test data.table.
#' @export

vl_motif_cl_enrich <- function(obj, 
                               cl_column)
{
  DT <- copy(obj)
  if(!cl_column %in% names(DT))
    stop("cl_column does not exist in provided obj")
  
  # Checks
  if(!all(c("motif", "motif_counts", "motif_name") %in% names(DT)))
    stop("Provided DT should contain c('motif', 'motif_counts', 'motif_name') columns. see ?vl_motif_counts() output!")
  
  # Rename cl column
  if(cl_column!="cl")
    names(DT)[names(DT)==cl_column] <- "cl"
  
  # Transform cl as factor and fill NAs with None
  if(!is.factor(DT$cl))
  {
    print("cl_column coerced as factor (levels in alphabetical order).")
    DT[, cl:= factor(cl, 
                     levels= sort(unique(DT$cl)))]
  }
  if(anyNA(DT$cl))
  {
    print("cl_column contains NAs that will be renamed as 'None' (factor)")
    DT[is.na(cl), cl:= factor("None")]
  }
  
  # Compute enrichment
  DT <- DT[, {
    .ccl <- data.table(V1= unique(cl))
    if(nrow(.ccl)>1)
      .ccl[, c("OR", "pval"):= {
        .f <- fisher.test(cl==V1, motif_counts>0)
        .f[c("estimate", "p.value")]
      }, V1]
  }, .(motif, motif_name)]
  DT[, padj:= p.adjust(pval, method = "fdr"), pval]
  DT[, log2OR:= log2(OR)]
  names(DT)[3] <- "cl"
  return(DT)
}
