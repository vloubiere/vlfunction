#' plot STRING interactions
#'
#' Network plotting funciton for vl_STRING_interaction output
#'
#' @param symbols vector of Dmel gene symbols
#' @param cex.vertices Scaling factor vertices
#' @param label.cex Scaling factor vertices labels
#' @param vertex.border.col Color of vertices' borders
#' 
#' @examples 
#' #USAGE
#' symbols <- c("Pc", "Psc", "E(z)", "dgt1", "Rcd1")
#' interactions <- vl_STRING_interaction(symbols)
#' vl_STRING_network(interactions)
#'
#' @return Network plot.
#' @export
#' 
vl_STRING_network <- function(obj,
                              cex.vertices= 1,
                              label.cex= 1,
                              vertex.border.col= "black")
{
  obj$vertices[, size:= size*cex.vertices]
  .g <- graph_from_data_frame(d = obj$edges, 
                              vertices = obj$vertices,
                              directed = F)
  plot(.g, 
       vertex.label.cex= label.cex, 
       vertex.frame.color= vertex.border.col)
}


