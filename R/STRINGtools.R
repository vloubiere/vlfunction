#' get STRING intrations DT
#'
#' Extract interactions form STRIN db (Dmel only for now)
#'
#' @param symbols vector of Dmel gene symbols
#' @param plot Should the graph be plotted?
#' @param score_cutoff If specified, only plot interactions with score>=cutoff
#' @param top_N If specified, only plot top N interactions
#' @param col For each gene symbol, the color of corresponding vertice (i.e, tomato for up, cornflowerblue for down)
#' @param size For each gene symbol, the size of corresponding vertice (i.e, absolute FC)
#' @param cex.label cex factor to be applied to correponding vertices labels
#' @examples 
#' #USAGE
#' symbols <- c("Pc", "Psc", "E(z)", "RpL10", "RpL11", "RpL12")
#' vl_STRING_interaction(symbols)
#' test <- vl_STRING_interaction(symbols, col= vl_palette_categ1(6))
#' plot(test)
#' plot(test, size= 10*1:6)
#' plot(test, col= "red")
#' @return An object that can be used with the plot.vl_STRING() method
#' @export
vl_STRING_interaction <- function(symbols,
                                  plot= T,
                                  score_cutoff= 900,
                                  top_N= NA,
                                  col= "tomato",
                                  size= 10,
                                  cex.label= 1)
{
  if(!identical(symbols, unique(symbols)))
    stop("symbols should be unique")
  
  # Get interaction
  res <- vl_Dmel_STRING_DB[protein1_symbol %in% symbols & protein2_symbol %in% symbols]
  
  # Vertices
  V <- data.table(name= symbols,
                  size, 
                  color= col,
                  cex.label)
  
  # Final obj
  obj <- list(symbols= symbols,
              interactions= res,
              vertices= V)
  class(obj) <- c("vl_STRING", "list")
  
  # Print info
  check <- 100-sum(symbols %in% c(vl_Dmel_STRING_DB$protein1_symbol, vl_Dmel_STRING_DB$protein2_symbol))/length(unique(symbols))*100
  print(paste0(check, "% of symbols have no correpondance in DB"))
  
  # PLOT
  if(plot)
    plot(obj,
         score_cutoff= score_cutoff,
         top_N= top_N,
         col= col,
         size= size,
         cex.label= cex.label)
  
  invisible(obj)
}

#' @describeIn vl_STRING_interaction Method to plot STRING interaction igraphs
#' @export
plot.vl_STRING <- function(obj,
                           score_cutoff= 900,
                           top_N= NA,
                           col,
                           size,
                           cex.label)
{
  list2env(obj, environment())
  
  # Checks and cutoffs
  if(!is.na(score_cutoff))
    interactions <- interactions[combined_score>=score_cutoff]
  if(!is.na(top_N))
    interactions <- interactions[order(combined_score)<=N_top]
  
  # Change vertices if new color, size or cex.label specified
  if(!missing(size))
    vertices$size <- size
  if(!missing(col))
    vertices$color <- col
  if(!missing(cex.label))
    vertices$cex.label <- cex.label
  
  # Remove vertices with no interactions
  vertices <- vertices[name %in% interactions[, c(protein1_symbol, protein2_symbol)]]
  
  # Make Edges object
  E <- interactions[, .(from= protein1_symbol,
                        to= protein2_symbol,
                        width= combined_score/999*3)]
  # igraph
  .g <- igraph::graph_from_data_frame(d = E, 
                                      vertices = vertices,
                                      directed = F)
  plot(.g, 
       vertex.label.cex= vertices$cex.label, 
       vertex.frame.color= NA)
}