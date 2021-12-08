#' plot STRING interactions
#'
#' Network plotting funciton for vl_STRING_interaction output
#'
#' @param symbols vector of Dmel gene symbols
#' @param cex.vertices Scaling factor vertices
#' @param cex.vertices.labels Scaling factor vertices labels
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
                              cex.vertices.labels= 1,
                              vertex.border.col= "black")
{
  if(length(cex.vertices)!=1)
    stop("lenght(cex.vertices) should be 1, applied to all vertices equally")
  if(length(cex.vertices.labels)!=1)
    stop("lenght(cex.vertices.labels) should be 1, applied to all vertices equally")
  
  .c <- data.table::copy(obj)
  .c$vertices[, size:= size*cex.vertices]
  .c$vertices[, cex.label:= cex.label*cex.vertices.labels]
  .g <- igraph::graph_from_data_frame(d = .c$edges, 
                                      vertices = .c$vertices,
                                      directed = F)
  plot(.g, 
       vertex.label.cex= .c$vertices$cex.label, 
       vertex.frame.color= vertex.border.col)
}


