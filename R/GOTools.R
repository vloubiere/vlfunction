#' GO analysis 
#'
#' This function compute GO enrichment in a group of genes
#'
#' @param FBgn_list A vector of FBgn IDs
#' @param go_object object containing GO data. see ?vl_Dmel_GO_FB2020_05
#' @param FBgn_universe Vector of FBgn IDs to which the analysis will be restricted. default= "all" (i.e. all the FBGns present in the table)
#' @param plot Plot result?
#' @param padj_cutoff cutoff for plotting
#' @param top_enrich Show only n top enriched motifs
#' @examples
#' RpL <- vl_genes_set[GO=="RpL_genes", FBgn]
#' par(mar= c(4,20,2,6))
#' vl_GO_enrich(RpL)
#' @export
vl_GO_enrich <- function(FBgn_vector,
                         go_object= vl_Dmel_GO_FB2020_05,
                         FBgn_universe= "all",
                         go_type= "all",
                         plot= T,
                         padj_cutoff= 0.05,
                         top_enrich= Inf)
{
  obj <- data.table::copy(go_object)
  FBgn_vector <- unique(FBgn_vector)
  # Filter DB with depending on Universe
  if(FBgn_universe!="all")
    obj <- obj[FBgn %in% FBgn_universe]
  # Keeps only GOs with at least one match
  obj <- obj[GO %chin% obj[FBgn_vector, on= "FBgn"][,GO]]
  # Counts genes per GO
  all_genes <- unique(obj$FBgn)
  check <- all_genes %chin% FBgn_vector
  obj[, c("OR", "pval"):= {
    tab <- table(all_genes %chin% FBgn,
                 check)
    fisher.test(tab, alternative = "greater")[c("estimate", "p.value")]
  }, GO]
  
  # padj...
  res <- obj[, .(GO, 
                 variable= name,
                 log2OR= log2(OR),
                 padj= p.adjust(pval, method = "fdr"))]
  res <- unique(res)
  setorderv(res, "log2OR")
  setattr(res, "class", c("vl_enr", "data.table", "data.frame"))

  #-----------------------#
  # Plot
  #-----------------------#
  if(plot)
    plot(res,
         padj_cutoff= padj_cutoff,
         top_enrich= top_enrich)

  invisible(res)  
}

#' GO analysis
#'
#' This function compares GO enrichments for several groups
#'
#' @param FBgn_list A nameds list of FBgn IDs (1 sublist/cluster)
#' @param go_object object containing GO data. see ?vl_fb_go_table_dm6_FB2020_05
#' @param all_FBgns Vector of FBgns to be used in the universe, typically all the genes tested with DESeq. Default= "all" FBgns present in go_object.
#' @param go_type The type of go to be considered. Default "all" means "biological_process" AND "cellular_component" AND "molecular_function". 
#' @param plot Should the result be plot using balloons plot?
#' @param padj_cutoff cutoff for ballons to be ploted
#' @param log2OR_cutoff cutoff for ballons to be ploted
#' @param N_top Select top enriched motifs/cluster
#' @param x_breaks Breaks used for ballon's sizes
#' @param color_breaks Color breaks used for coloring
#' @param cex.balloons Expansion factor for balloons
#' @param col Vector of colors used for coloring
#' @param main Title. Default= NA
#' @param auto_margins Use auto margins? Default= T
#' @examples
#' FBgn_list <- split(vl_genes_set$FBgn, vl_genes_set$GO)
#' par(mar= c(5,25,2,5), las= 1)
#' vl_GO_clusters(FBgn_list, auto_margins = F, cex.balloons = 0.3)
#' @export
vl_GO_clusters <- function(FBgn_list,
                           go_object= vl_Dmel_GO_FB2020_05,
                           FBgn_universe= "all",
                           go_type= "all",
                           plot= T,
                           padj_cutoff= 0.00001,
                           log2OR_cutoff= 0,
                           N_top= Inf,
                           x_breaks,
                           color_breaks,
                           cex.balloons= 1,
                           col= c("cornflowerblue", "lightgrey", "tomato"),
                           main= NA,
                           auto_margins = T)
{
  obj <- data.table::copy(go_object)
  # Filter DB with depending on Universe
  if(FBgn_universe!="all")
    obj <- obj[FBgn %in% FBgn_universe]
  if(is.null(names(FBgn_list)))
    names(FBgn_list) <- paste0("cl", seq(FBgn_list))
  
  # Compute enrichment per cluster
  obj <- lapply(FBgn_list, function(x) {
    vl_GO_enrich(x,
                 go_object,
                 FBgn_universe,
                 go_type,
                 plot= F)
  })
  names(obj) <- names(FBgn_list)
  
  # Result
  res <- rbindlist(obj, idcol = "cl")
  class(res) <- c("vl_enr_cl", "data.table", "data.frame")
  
  #---------------------------------------#
  # Plot
  #---------------------------------------#
  if(plot)
    plot(res,
         padj_cutoff= padj_cutoff,
         log2OR_cutoff= log2OR_cutoff,
         N_top= N_top,
         x_breaks= x_breaks,
         color_breaks= color_breaks,
         col= col,
         main= main,
         cex.balloons= cex.balloons,
         auto_margins = auto_margins)
  
  invisible(res)
}