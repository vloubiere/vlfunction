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
#' @param x_breaks Breaks used for balloons sizes
#' @param color_breaks Balloons color breaks
#' @param cex.balloons Balloons size expansion factor
#' @param col Color vector
#' @param main Title
#' @param plot_empty_clusters Should empty clusters be plotted? Default to TRUE
#' @param order Value to be used for ordering before selecting top enriched. Possible values are "padj", "log2OR". defaut= "padj"
#' @examples
#' RpL <- vl_genes_set[GO=="RpL_genes", FBgn]
#' Hox <- vl_genes_set[GO=="HOX_genes", FBgn]
#' par(mar= c(4,20,2,6), las= 1)
#' vl_GO_enrich(RpL, species= "Dm", plot= T)
#' DT <- vl_GO_enrich(Hox, species= "Dm")
#' plot(DT, padj_cutoff= 1e-10)
#' DT <- vl_GO_enrich(list(RpL, Hox), species= "Dm", cex.balloons = 0.3, padj_cutoff= 1e-5, top_enrich = 20, plot= T)
#' plot(DT, top_enrich= 5)
#' @export
vl_GO_enrich <- function(geneIDs,
                         geneUniverse_IDs= NULL,
                         species,
                         plot= F,
                         padj_cutoff= 1e-5,
                         top_enrich= NA,
                         order= "padj",
                         x_breaks,
                         color_breaks,
                         cex.balloons= 1,
                         col= c("blue", "red"),
                         main= NA,
                         plot_empty_clusters= T)
{
  # Checks
  if(is.character(geneIDs))
    geneIDs <- list(set= geneIDs)
  if(is.null(names(geneIDs)))
    names(geneIDs) <- seq(geneIDs)
  if(anyDuplicated(names(geneIDs)))
    stop("names(geneIDs) should be unique!")
  # Make unique
  geneIDs <- lapply(geneIDs, unique)
  
  # Genome
  db <- switch(species,
               "Dm"= org.Dm.eg.db::org.Dm.eg.db,
               "Mm"= org.Mm.eg.db::org.Mm.eg.db)
  keyType <- switch(species, 
                    "Dm"= "FLYBASE",
                    "Mm"= "ENSEMBL")
  
  #-----------------------------------#
  # Extract sets and universe GOs
  #-----------------------------------#
  # Sets
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
  if(!is.null(geneUniverse_IDs) && any(!set$ID %in% geneUniverse_IDs))
    stop("Some geneIDs are not included in the Universe!")
  # Universe
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
  DT <- CJ(variable= unique(set$GO), 
           cl= factor(names(geneIDs), names(geneIDs)))
  # Overlaps set
  setkeyv(set, "GO")
  DT[, set_hit:= sum(geneIDs[[cl]] %chin% set[.BY[2], ID]), .(cl, variable)] # Genes from cluster found in GO
  DT[, set_total:= sum(geneIDs[[cl]] %chin% set$ID), cl] # Total genes in cluster (existing in GO database)
  # Overlaps universe
  setkeyv(uni, "GO")
  DT[, ctl_hit:= uni[.BY, .N], variable]# No cluster there= simpler!
  DT[, ctl_total:= length(unique(uni$ID))]# Total universe
  # Fisher test
  DT[, c("OR", "pval"):= {
    fisher.test(matrix(unlist(.BY), byrow= T, ncol= 2))[c("estimate", "p.value")]
  }, .(set_hit, set_total, ctl_hit, ctl_total)]
  
  ###############################
  # Add GO description and clean
  ###############################
  terms <- AnnotationDbi::select(GO.db::GO.db,
                                 keys = as.character(unique(DT$variable)),
                                 keytype= "GOID", 
                                 columns= "TERM")
  terms <- as.data.table(terms)
  DT[terms, name:= i.TERM, on= "variable==GOID"]
  DT[, log2OR:= log2(OR)]
  DT[, padj:= p.adjust(pval, method = "fdr"), cl]
  DT$OR <- DT$pval <- NULL
  DT <- na.omit(DT)
  setorderv(DT, c("cl", "log2OR"))

  ###############################
  # Plot and export
  ###############################
  if(length(levels(DT$cl))>1)
  {
    setattr(DT, "class", c("vl_enr_cl", "data.table", "data.frame"))
    if(plot)
      DT <- plot.vl_enr_cl(obj= DT,
                           padj_cutoff= padj_cutoff,
                           top_enrich= top_enrich, 
                           order= order,
                           x_breaks= x_breaks,
                           color_breaks= color_breaks,
                           cex.balloons= cex.balloons,
                           col= col,
                           main= main, 
                           plot_empty_clusters = plot_empty_clusters)
  }else{
    setattr(DT, "class", c("vl_enr", "data.table", "data.frame"))
    if(plot)
    {
      plot.vl_enr(obj= DT,
                  padj_cutoff= padj_cutoff,
                  top_enrich= top_enrich,
                  order= order,
                  xlab = "Odd Ratio (log2)",
                  col = col)
      if(!is.na(main))
        title(main= main)
    }
  }

  invisible(DT)  
}
