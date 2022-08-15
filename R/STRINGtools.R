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
                                  score_cutoff= 400,
                                  top_N= Inf,
                                  col= "tomato",
                                  size= 10,
                                  cex.label= 1)
{
  string_db <- switch(species,
                      "Dm" = STRINGdb::STRINGdb$new(version = "10", species = 7227, score_threshold = 0, input_directory = ""),
                      "Mm" = STRINGdb::STRINGdb$new(version = "10", species = 1090, score_threshold = 0, input_directory = ""))
  # Vertices
  V <- data.table(name= symbols,
                  size, 
                  color= col,
                  cex.label)
  
  # Get IDs (STRING put them in CAPS!)
  IDs <- as.data.table(string_db$map(V, "name", removeUnmappedRows = TRUE))
  IDs[V[, .(name, caps= toupper(name))], name:= i.name, on= "name==caps"]
  
  # Full graphs
  E <- string_db$get_interactions(IDs$STRING_id)
  E <- as.data.table(E)
  E[IDs, from:= i.name, on= "from==STRING_id"]
  E[IDs, to:= name, on= "to==STRING_id"]
  E <- E[, .(width= max(combined_score)), .(from, to)]
  
  # Keep only V and E interacting within subgroup
  V <- V[name %in% E[, c(from, to)]]
  E <- E[from %in% V$name & to %in% V$name]
  setorderv(E, "width", -1)
  
  # igraph
  obj <- list(V= V, E= E)
  class(obj) <- c("vl_STRING", "list")
  
  # PLOT
  if(plot)
    plot(obj= obj, 
         score_cutoff= score_cutoff,
         top_N= top_N)
  
  invisible(obj)
}

#' @describeIn vl_STRING_interaction Method to plot STRING interaction igraphs
#' @export
plot.vl_STRING <- function(obj,
                           score_cutoff= 400,
                           top_N= Inf)
{
  E <- obj$E
  V <- obj$V
  
  # Checks and cutoffs
  E <- E[width>=score_cutoff]
  E <- E[seq(nrow(E))<=top_N]
  E[, width:= width/999*3]
  
  # Keep only V and E interacting within subgroup
  V <- V[name %in% E[, c(from, to)]]
  E <- E[from %in% V$name & to %in% V$name,]
  
  # igraph
  .g <- igraph::graph_from_data_frame(d = E, 
                                      vertices = V,
                                      directed = F)
  plot(.g, 
       vertex.label.cex= V$cex.label,
       vertex.frame.color= NA)
}
