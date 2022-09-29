#' GO analysis 
#'
#' This function compute GO enrichment in a group of genes
#'
#' @param geneIDs A vector (barplot) of a list (clusters) of gene IDs
#' @param geneUniverse_IDs Vector of FBgn IDs to which the analysis will be restricted. If NULL, uses all genes.
#' @param species Chose between "Dm" and "Mm" 
#' @param plot Plot result?
#' @param padj_cutoff cutoff for plotting
#' @param top_enrich Show only n top enriched motifs
#' @param breaks Color breaks to be used. Defaults to range of filtered padj.
#' @examples
#' RpL <- vl_genes_set[GO=="RpL_genes", FBgn]
#' Hox <- vl_genes_set[GO=="HOX_genes", FBgn]
#' par(mar= c(4,20,2,6))
#' vl_GO_enrich(RpL, species= "Dm")
#' vl_GO_enrich(Hox, species= "Dm")
#' vl_GO_enrich(list(RpL, Hox), species= "Dm", cex.balloons = 0.3, padj_cutoff= 1e-5, top_enrich = 20)
#' @export
vl_GO_enrich <- function(geneIDs,
                         geneUniverse_IDs= NULL,
                         species,
                         plot= T,
                         padj_cutoff= 0.05,
                         top_enrich= Inf,
                         breaks= NULL)
{
  db <- switch(species,
               "Dm"= org.Dm.eg.db::org.Dm.eg.db,
               "Mm"= org.Mm.eg.db::org.Mm.eg.db)
  keyType <- switch(species, 
                    "Dm"= "FLYBASE",
                    "Mm"= "ENSEMBL")
  ###############################
  # Extract sets and universe GOs
  ###############################
  set <- AnnotationDbi::select(x= db,
                              keys = unique(unlist(geneIDs)),
                              keytype= keyType, 
                              columns= "GOALL")
  set <- data.table::as.data.table(set)
  setnames(set, 
           c(keyType, "GOALL"), 
           c("ID", "GO"))
  set <- unique(set[, .(ID, GO)]) # GOs IDs are reported for each evidence type
  set <- na.omit(set)
  uni <- AnnotationDbi::select(x= db,
                               keys = unique(set$GO), # Only GOs from test set are relevant!
                               keytype= "GOALL", 
                               columns= keyType)
  uni <- data.table::as.data.table(uni)
  setnames(uni, 
           c(keyType, "GOALL"),
           c("ID", "GO"))
  uni <- unique(uni[, .(ID, GO)]) # GOs IDs are reported for each evidence type
  if(!is.null(geneUniverse_IDs))
    uni <- uni[ID %chin% geneUniverse_IDs]
  ###############################
  # Format objects and Compute enrichments
  ###############################
  if(is.character(geneIDs))
    geneIDs <- list(set= geneIDs)
  DT <- rbindlist(lapply(geneIDs, function(x) SJ(ID= x)), idcol = "cl")
  DT <- unique(DT)
  DT[, cl:= factor(cl)]
  setkeyv(set, "ID")
  setkeyv(uni, "GO")
  res <- DT[, {
    # Extract genes from current cl
    genes <- ID
    genes <- set[genes]
    genes[, set:= T]
    # Extract universe genes with matching GOs
    bg <- unique(genes$GO)
    bg <- uni[bg]
    bg[, set:= F]
    # rbind, dcast, melt to get counts/gene/GO
    .c <- rbind(genes[, .(ID, GO, set)], 
                bg[, .(ID, GO, set)])
    .c <- dcast(.c, 
                set+ID~GO, 
                fun.aggregate = length)
    .c <- melt(.c, id.vars = c("set", "ID"), variable.name = "GO")
    # Fisher test
    .c[, {
      tab <- table(factor(set, c(T,F)),
                   factor(value>0, c(T,F)))
      fisher.test(tab, alternative = "greater")[c("estimate", "p.value")]
    }, GO]
  }, cl]
  ###############################
  # Add GO description and clean
  ###############################
  terms <- AnnotationDbi::select(GO.db::GO.db,
                                 keys = as.character(unique(res$GO)),
                                 keytype= "GOID", 
                                 columns= "TERM")
  terms <- as.data.table(terms)
  res[terms, variable:= i.TERM, on= "GO==GOID"]
  res[, log2OR:= log2(estimate)]
  res[, padj:= p.adjust(p.value, method = "fdr"), cl]
  res$estimate <- NULL
  res$p.value <- NULL
  res <- na.omit(res)
  setorderv(res, c("cl", "log2OR"))

  ###############################
  # Plot and export
  ###############################
  if(length(unique(res$cl))>1)
  {
    setattr(res, "class", c("vl_enr_cl", "data.table", "data.frame"))
    if(plot)
      plot(res,
           padj_cutoff= padj_cutoff,
           top_enrich= top_enrich)
  }else{
    setattr(res, "class", c("vl_enr", "data.table", "data.frame"))
    if(plot)
      plot(res,
           padj_cutoff= padj_cutoff,
           top_enrich= top_enrich)
  }

  invisible(res)  
}
