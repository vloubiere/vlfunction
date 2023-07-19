#' Build STRING db
#' @param species either "Dm" or "Mm"
#' @param network_type The type of interactions to be used. Can be one of "full" (full functional, default) or "physical" if only physical interactions are to be considered
#' @param version databse version. default= "10"
#' 
#' @export
vl_STRING_getDB <- function(species,
                            network_type= "full",
                            version= "10")
{
  if(!network_type %in% c("full", "physical"))
    stop("Network_type has to be one of 'full' (full functional annot) or 'physical' (physical interactions)")
  
  STRINGdb::STRINGdb$new(version = version, 
                         species = switch(species, "Dm"= 7227, "Mm"= 10090),
                         score_threshold = 0,
                         network_type = network_type,
                         input_directory = "")
}

#' get STRING interactions Object
#'
#' Extract interactions form STRINGdb
#'
#' @param STRINGdb A STRINGdb object. See ?vl_STRING_getDB
#' @param symbols vector of Dmel gene symbols
#' @param score_cutoff If specified, only plot interactions with score>=cutoff
#' @param top_N If specified, only plot top N interactions
#' @param remove_non_connected Remove proteins that have 0 connections. Default= T
#' @param plot Should the graph be plotted?
#' @param size Vertex size. Default= 15
#' @param size2 Vertex size2 for certain shapes. Default= 15
#' @param color Vertex color. Default= "tomato"
#' @param frame.color Frame color. Default= "black"
#' @param frame.width Frame width. Default= 1
#' @param shape Shape of the vertex. Default= "circle"
#' @param label Label of the vertex. Default= symbols
#' @param label.family Default= "serif"
#' @param label.font Default= 1
#' @param label.cex Default= 1
#' @param label.dist The distance of the label from the center of the vertex. If it is 0 then the label is centered on the vertex. If it is 1 then the label is displayed beside the vertex. Default= 0
#' @param label.degree It defines the position of the vertex labels, relative to the center of the vertices. It is interpreted as an angle in radian. Default= -pi/4.
#' @param label.color Default= "black"
#' @examples
#' # Build database
#' db <- vl_STRING_getDB(species= "Dm", network_type = "full", version= "11")
#' 
#' # Extract interactions
#' symbols <- c("Pc", "Psc", "E(z)", "RpL10", "RpL11", "RpL12")
#' vl_STRING_interaction(STRINGdb = db,
#'                       symbols = symbols)
#' .i <- vl_STRING_interaction(STRINGdb = db,
#'                             symbols = symbols,
#'                             plot= F)
#' plot(.i)
#' 
#' #Modify ploting parameters
#' print(paste0("param= ", paste0(colnames(igraph::as_data_frame(.i, what= "vertices")), collapse = ", ")))
#' igraph::V(.i)$color <- "cornflowerblue"
#' igraph::V(.i)$label.cex <- c(1,1,2,2,3,3)
#' plot(.i)
#' 
#' 
#' .i$V$color <- "cornflowerblue"
#' .i$V$label.cex <- seq(1, 3, length.out= 6)
#' .i$V$size <- seq(10, 60, 10)
#' plot(.i)
#' 
#' @return An igraph network. See ?plot.igraph()
#' @export
vl_STRING_interaction <- function(STRINGdb,
                                  symbols,
                                  score_cutoff= 400,
                                  top_N= Inf,
                                  remove_non_connected= T,
                                  plot= T,
                                  size= 15,
                                  size2= 15,
                                  color= "tomato",
                                  frame.color= "black",
                                  frame.width= 1,
                                  shape= "circle",
                                  label= symbols,
                                  label.family= "serif",
                                  label.font= 1,
                                  label.cex= 1,
                                  label.dist= 0,
                                  label.degree= -pi/4,
                                  label.color= "black")
{
  # Vertices
  V <- data.table(symbol= symbols,
                  size= size,
                  size2= size2,
                  color= color,
                  frame.color= frame.color,
                  frame.width= frame.width,
                  shape= shape,
                  label= label,
                  label.family= label.family,
                  label.font= label.font,
                  label.cex= label.cex,
                  label.dist= label.dist,
                  label.degree= label.degree,
                  label.color= label.color)
  V <- unique(V)
  
  # Get IDs (STRING put them in CAPS!)
  IDs <- as.data.table(STRINGdb$map(V, "symbol", removeUnmappedRows = TRUE))
  IDs[V[, .(symbol, caps= toupper(symbol))], symbol:= symbol, on= "symbol==caps"]
  
  # Full graphs
  E <- STRINGdb$get_interactions(IDs$STRING_id)
  E <- as.data.table(E)
  E[IDs, from:= i.symbol, on= "from==STRING_id"]
  E[IDs, to:= symbol, on= "to==STRING_id"]
  E <- E[, .(width= max(combined_score)), .(from, to)]
  
  # Checks and cutoffs
  E <- E[width>=score_cutoff]
  setorderv(E, "width", -1)
  E <- E[seq(nrow(E))<=top_N]
  E[, width:= width/999*3]
  
  # Keep only V and E interacting within subgroup
  if(remove_non_connected)
  {
    V <- V[symbol %in% E[, c(from, to)]]
    E <- E[from %in% V$symbol & to %in% V$symbol]
  }
  
  # Make igraph object
  .i <- igraph::graph_from_data_frame(d = E, 
                                      vertices = V,
                                      directed = F)
  
  # PLOT
  if(plot)
    plot(.i)
  
  # Return
  invisible(.i)
}

