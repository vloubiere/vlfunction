#' get STRING intrations DT
#'
#' Extract interactions form STRIN db (Dmel only for now)
#'
#' @param symbols vector of Dmel gene symbols
#' @param species either "Dm" or "Mm"
#' @param plot Should the graph be plotted?
#' @param score_cutoff If specified, only plot interactions with score>=cutoff
#' @param top_N If specified, only plot top N interactions
#' @param col For each gene symbol, the color of corresponding vertice (i.e, tomato for up, cornflowerblue for down)
#' @param size For each gene symbol, the size of corresponding vertice (i.e, absolute FC)
#' @param cex.label cex factor to be applied to correponding vertices labels
#' @examples 
#' #USAGE
#' symbols <- c("Pc", "Psc", "E(z)", "RpL10", "RpL11", "RpL12")
#' net <- vl_STRING_interaction(symbols, species= "Dm", col= vl_palette_few_categ(6))
#' E <- igraph::as_data_frame(net, what= "edges") 
#' V <- igraph::as_data_frame(net, what= "vertices") 
#' V$color <- grey.colors(6)
#' .g <- igraph::graph_from_data_frame(d = E, 
#'                                     vertices = V,
#'                                     directed = F)
#' plot(.g)
#' @return An igraph network
#' @export
vl_STRING_interaction <- function(symbols,
                                  species,
                                  plot= T,
                                  score_cutoff= 900,
                                  top_N= NA,
                                  col= "tomato",
                                  size= 10,
                                  cex.label= 1)
{
  string_db <- switch(species,
                      "Dm" = STRINGdb$new(version = "10", species = 7227, score_threshold = score_cutoff, input_directory = ""),
                      "Mm" = STRINGdb$new(version = "10", species = 1090, score_threshold = score_cutoff, input_directory = ""))
  V <- data.table(symbol= symbols,
                   size, 
                   color= col,
                   cex.label)
  V <- as.data.table(string_db$map(V, "symbol", removeUnmappedRows = TRUE))
  net <- string_db$get_subnetwork(V$STRING_id)
  
  E <- igraph::as_data_frame(net, what= "edges")
  E <- as.data.table(E)
  E <- E[, .(from, to, width= combined_score/999*3)] 
  E[V, from:= symbol, on= "from==STRING_id"]
  E[V, to:= symbol, on= "to==STRING_id"]
  
  # igraph
  .g <- igraph::graph_from_data_frame(d = E, 
                                      vertices = V,
                                      directed = F)
  
  # PLOT
  if(plot)
    plot(.g, 
         vertex.label.cex= vertices$cex.label, 
         vertex.frame.color= NA)
  
  
  invisible(.g)
}
