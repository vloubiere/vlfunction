#' get STRING intrations DT
#'
#' Extract interactions form STRIN db (Dmel only for now)
#'
#' @param symbols vector of Dmel gene symbols
#' @param size For each gene symbol, the size of corresponding vertice (i.e, absolute FC)
#' @param col For each gene symbol, the color of corresponding vertice (i.e, tomato for up, cornflowerblue for down)
#' @param score_cutoff Interaction score cutoff. Default is 200, max value in DB is 999. Also used to compute edge widths!
#' @param db_path Path to an existing db, or filename where to create it. See examples
#' @examples 
#' #USAGE
#' symbols <- c("Pc", "Psc", "E(z)", "dgt1", "Rcd1")
#' interactions <- vl_STRING_interaction(symbols)
#' vl_STRING_network(interactions)
#' 
#' 
#' #------------#
#' # Code used to create DB
#' #------------#
#' # Download STRING db v11.5 fro Dmel
#' link <- "https://stringdb-static.org/download/protein.links.detailed.v11.5/7227.protein.links.detailed.v11.5.txt.gz"
#' tmp <- tempfile(fileext = ".txt.gz")
#' download.file(link, tmp)
#' DB <- fread(tmp)
#' # Download correspondance Dmel symbols
#' link <- "https://stringdb-static.org/download/protein.info.v11.5/7227.protein.info.v11.5.txt.gz"
#' download.file(link, tmp)
#' sym <- fread(tmp)
#' DB[sym, protein1_symbol:= i.preferred_name, on= "protein1==`#string_protein_id`"]
#' DB[sym, protein2_symbol:= i.preferred_name, on= "protein2==`#string_protein_id`"]
#' # Collapse unique interactions
#' DB[, merge:= paste0(sort(c(protein1_symbol, protein2_symbol)), collapse= ""), .(protein1_symbol, protein2_symbol)]
#' DB <- DB[, .SD[1], merge]
#' # SAVE
#' saveRDS(DB[, .(protein1_symbol, protein2_symbol, combined_score)], 
#' file = db_path)
#' 
#' @return An object that can be used with the vl_STRING_network() function.
#' @export

vl_STRING_interaction <- function(symbols= NULL, 
                                  size= NULL,
                                  col= NULL,
                                  score_cutoff= 200,
                                  top_N= NA,
                                  db_path= "/mnt/d/_R_data/genomes/dm6/STRING/dmel_7227_v11.5_vl_db.rds")
{
  # Checks
  if(is.null(size))
    size <- rep(20, length(symbols))
  if(is.null(col))
    col <- rep("tomato", length(symbols))
  # Import DB
  if(is.na(db_path))
    stop("db_path should either be a path to a STRING_db object (vl) or a valid filename where the db will be created")
  if(!file.exists(db_path))
  {
    check <- readline("db_path file does not exist. Should the db be created in place (~500M)? y/n ")
    if(check=="y")
    {
      if(!grepl(".rds$", db_path))
        stop("db_path should point to a valid/accesssible .rds filename")
      # Download STRING db v11.5 fro Dmel
      link <- "https://stringdb-static.org/download/protein.links.detailed.v11.5/7227.protein.links.detailed.v11.5.txt.gz"
      tmp <- tempfile(fileext = ".txt.gz")
      download.file(link, tmp)
      DB <- fread(tmp)
      # Download correspondance Dmel symbols
      link <- "https://stringdb-static.org/download/protein.info.v11.5/7227.protein.info.v11.5.txt.gz"
      download.file(link, tmp)
      sym <- fread(tmp)
      DB[sym, protein1_symbol:= i.preferred_name, on= "protein1==`#string_protein_id`"]
      DB[sym, protein2_symbol:= i.preferred_name, on= "protein2==`#string_protein_id`"]
      # Collapse unique interactions
      DB[, merge:= paste0(sort(c(protein1_symbol, protein2_symbol)), collapse= ""), .(protein1_symbol, protein2_symbol)]
      DB <- DB[, .SD[1], merge]
      # SAVE
      saveRDS(DB[, .(protein1_symbol, protein2_symbol, combined_score)], 
              file = db_path)
    }
  }else
    DB <- readRDS(db_path)
  
  # Extract interactions
  sub <- DB[protein1_symbol %in% symbols & 
              protein2_symbol %in% symbols]
  # Print info
  check <- 100-length(which(symbols %in% c(DB$protein1, DB$protein2)))/length(unique(symbols))*100
  print(paste0(check, "% of symbols have not correpondance in DB"))
  # Score cutoff
  sub <- sub[combined_score>=score_cutoff]
  # Make Vertices object
  V <- data.table(name= symbols,
                  size, 
                  color= col)
  V <- unique(V)
  # Make Edges object
  E <- V[, sub[.BY, .(to= protein2_symbol, 
                      width= combined_score/999*5), on= "protein1_symbol==from"], .(from= name)]
  E <- na.omit(E)
  setorderv(E, "width", order = -1)
  if(!is.na(top_N) & nrow(E)>top_N)
    E <- E[1:top_N]
  # Keep only vertices which are represented
  V <- V[name %in% E$from | name %in% E$to]
  # RETURN
  return(list(vertices= V, 
              edges= E))
}
