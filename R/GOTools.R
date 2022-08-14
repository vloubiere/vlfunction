#' GO analysis 
#'
#' This function compute GO enrichment in a group of genes
#'
#' @param geneSet_IDs A vector of FBgn IDs
#' @param geneUniverse_IDs Vector of FBgn IDs to which the analysis will be restricted. If NULL, uses all genes.
#' @param species Chose between "Dm" and "Mm" 
#' @param plot Plot result?
#' @param padj_cutoff cutoff for plotting
#' @param top_enrich Show only n top enriched motifs
#' @examples
#' RpL <- vl_genes_set[GO=="RpL_genes", FBgn]
#' Hox <- vl_genes_set[GO=="HOX_genes", FBgn]
#' par(mar= c(4,20,2,6))
#' vl_GO_enrich(RpL, species= "Dm")
#' vl_GO_enrich(Hox, species= "Dm")
#' vl_GO_enrich(RpL, Hox, species= "Dm")
#' @export
vl_GO_enrich <- function(geneSet_IDs,
                         geneUniverse_IDs= NULL,
                         species,
                         plot= T,
                         padj_cutoff= 0.05,
                         top_enrich= Inf)
{
  db <- switch(species,
               "Dm"= org.Dm.eg.db::org.Dm.eg.db,
               "Mm"= org.Mm.eg.db::org.Mm.eg.db)
  keyType <- switch(species, 
                    "Dm"= "FLYBASE",
                    "Mm"= "ENSEMBL")
  # Extract geneSet GOs
  set <- AnnotationDbi::select(x= db,
                               keys = unique(c(geneSet_IDs)),
                               keytype= keyType, 
                               columns= "GOALL")
  set <- data.table::as.data.table(set)
  setnames(set, 
           c(keyType, "GOALL"), 
           c("ID", "GO"))
  set <- unique(set[, .(ID, GO)])
  # Extract GOs genes and restrict to universe
  uni <- AnnotationDbi::select(x= db,
                               keys = unique(set$GO),
                               keytype= "GOALL", 
                               columns= keyType)
  uni <- data.table::as.data.table(uni)
  setnames(uni, 
           c(keyType, "GOALL"),
           c("ID", "GO"))
  if(!is.null(geneUniverse_IDs))
    uni <- uni[set$ID, on= "ID"]
  # Format and compute
  DT <- rbind(SJ(set[, .(GO, ID)], set= T),
              SJ(uni[, .(GO, ID)], set= F))
  DT <- dcast(DT, set+ID~GO, fun.aggregate = length)
  DT <- melt(DT, id.vars = c("set", "ID"), variable.name = "GO")
  res <- DT[, {
    tab <- table(factor(set, c(T,F)),
                 factor(value>0, c(T,F)))
    fisher.test(tab, alternative = "greater")[c("estimate", "p.value")]
  }, GO]
  res[, variable:= AnnotationDbi::select(GO.db::GO.db,
                                         keys = as.character(GO),
                                         keytype= "GOID", 
                                         columns= "TERM")$TERM]
  # padj and export
  res[, c("log2OR", "padj"):= .(log2(estimate), p.adjust(p.value, method = "fdr"))]
  res$estimate <- NULL
  res$p.value <- NULL
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
#' @param FBgn_list A named list of FBgn IDs (1 sublist/cluster)
#' @param geneUniverse_IDs Vector of FBgns to be used in the universe, typically all the genes tested with DESeq. If NULL, uses all genes.
#' @param species Chose between "Dm" and "Mm" 
#' @param plot Should the result be plot using balloons plot?
#' @param padj_cutoff cutoff for ballons to be ploted
#' @param log2OR_cutoff cutoff for ballons to be ploted
#' @param top_enrich Select top enriched motifs/cluster
#' @param x_breaks Breaks used for ballon's sizes
#' @param color_breaks Color breaks used for coloring
#' @param cex.balloons Expansion factor for balloons
#' @param col Vector of colors used for coloring
#' @param main Title. Default= NA
#' @param auto_margins Use auto margins? Default= T
#' @examples
#' FBgn_list <- split(vl_genes_set$FBgn, vl_genes_set$GO)
#' par(mar= c(5,25,2,5), las= 1)
#' vl_GO_clusters(FBgn_list, species= "Dm", auto_margins = F, cex.balloons = 0.3, padj_cutoff= 1e-5, top_enrich = 20)
#' @export
vl_GO_clusters <- function(FBgn_list,
                           geneUniverse_IDs= NULL,
                           species,
                           plot= T,
                           padj_cutoff= 0.05,
                           log2OR_cutoff= 0,
                           top_enrich= Inf,
                           x_breaks,
                           color_breaks,
                           cex.balloons= 1,
                           col= c("cornflowerblue", "lightgrey", "tomato"),
                           main= NA,
                           auto_margins = T)
{
  # Compute enrichment per cluster
  obj <- lapply(FBgn_list, function(x) {
    vl_GO_enrich(x,
                 geneUniverse_IDs,
                 species= species,
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
         top_enrich= top_enrich,
         x_breaks= x_breaks,
         color_breaks= color_breaks,
         col= col,
         main= main,
         cex.balloons= cex.balloons,
         auto_margins = auto_margins)
  
  invisible(res)
}
