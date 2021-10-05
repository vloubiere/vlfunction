#' Compute motif enrichment
#'
#' Compute motif enrichment for the cluster in cl_columns, using all the lines as background
#'
#' @param obj A motif count object similar to the output of ?vl_motif_counts() (colnames:'motif', 'motif_counts', 'motif_name') + one cluster column! 
#' @param cl_column name of the cluster column
#' @param bg cluster IDs to be used as background. Default unique(obj[[cl_column]])
#' @param comp_expr Expression used for contingency table. Default "motif_count>0"
#' 
#' @examples 
#' 
#' @return Fisher test data.table.
#' @export

vl_motif_cl_enrich <- function(obj, 
                               cl_column,
                               bg= unique(obj[[cl_column]]),
                               comp_expr= "motif_counts>0")
{
  DT <- copy(obj)
  # Checks
  if(!cl_column %in% names(DT))
    stop("cl_column does not exist in provided obj") else if(cl_column != "cl")
      names(DT)[names(DT)==cl_column] <- "cl"
  if(class(DT$cl) != class(bg))
    stop("cluster column and bg arg should share the same class!")
  if(!all(c("motif", "motif_counts", "motif_name") %in% names(DT)))
    stop("Provided DT should contain c('motif', 'motif_counts', 'motif_name') columns. see ?vl_motif_counts() output!")
  if(anyNA(bg))
    stop("bg contains NAs. cl_column should not contain NAs OR they should not be used for comparisons!")
  DT <- DT[, .(motif, motif_name, motif_counts, cl)]
  setkeyv(DT, c("motif", "cl"))
  N_mot <- length(unique(DT$motif))
    
  # Format result table
  cmb <- DT[, {
    data.table(V1= na.omit(unique(cl))) # make DT containing clusters to test
  }, .(mot= motif, motif_name)]
  
  # For each motif/cluster combination, compute association using fisher
  cmb[, c("OR", "pval"):= {
    # Extract motif from DT, restrict to regions from tested cluster OR bg, cast contingency table
    .t <- table(DT[.(mot, c(V1, bg)), .(cl==V1, 
                                        motif_counts>0)])
    # If motif present in both tested cl and bg, do fisher test
    if(identical(dim(.t), as.integer(c(2,2))))
      res <- as.data.table(as.list(fisher.test(.t)[c("estimate", "p.value")])) else
        data.table(numeric(), numeric())
  }, .(mot, V1)] 
  cmb <- na.omit(cmb)
  cmb[, padj:= p.adjust(pval, method = "fdr"), pval]
  cmb[, log2OR:= log2(OR)]
  names(cmb)[c(1,3)] <- c("motif", cl_column)
  return(cmb)
}
