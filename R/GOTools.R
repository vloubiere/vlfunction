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
  # Select
  counts <- go_object$FBgn_GO_counts_matrix
  if(FBgn_universe != "all")
    counts <- counts[FBgn %in% FBgn_universe]
  names <- go_object$GO_names
  if(go_type != "all")
  {
    names <- names[type %in% go_type]
    cols <- c("FBgn", names$GO)
    counts <- counts[, ..cols]
  }

  # Select GOs for which at least one of the genes is represented 
  sel <- counts[FBgn_vector, apply(.SD, 2, function(x) any(x>0)), on= "FBgn", .SDcols= patterns("GO:")]
  sel <- c("FBgn", names(sel)[sel])
  counts <- counts[, ..sel]
  # Melt
  .m <- melt(counts, id.vars = "FBgn")
  .m[FBgn_vector, cl:= 1, on= "FBgn"]
  .m[is.na(cl), cl:= 0]
  # Compute enrichment
  res <- .m[, {
    tab <- table(value>0, cl>0)
    if(identical(dim(tab), c(2L,2L)))
    {
      .f <- fisher.test(tab, 
                        alternative = "greater")
      .(OR= .f$estimate,
        pval= .f$p.value)
    }else
      .(OR= as.numeric(NA),
        pval= as.numeric(NA))
  }, variable]
  
  #-----------------------#
  # Plot
  #-----------------------#
  # padj...
  res[, padj:= p.adjust(pval, method = "fdr"), pval]
  res[, log2OR:= log2(OR)]
  res[names, variable:= i.name, on= "variable==GO"]
  res <- res[, .(variable, log2OR, padj)]
  setorderv(res, "log2OR")
  setattr(res, "class", c("vl_enr", "data.table", "data.frame"))
  
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
  counts <- go_object$FBgn_GO_counts_matrix
  if(FBgn_universe != "all")
    counts <- counts[FBgn %in% FBgn_universe]
  names <- go_object$GO_names
  if(go_type != "all")
  {
    names <- names[type %in% go_type]
    cols <- c("FBgn", names$GO)
    counts <- counts[, ..cols]
  }
  
  # Select GOs for which at least one of the genes is represented 
  sel <- counts[unique(unlist(FBgn_list)), apply(.SD, 2, function(x) any(x>0)), on= "FBgn", .SDcols= patterns("GO:")]
  sel <- c("FBgn", names(sel)[sel])
  counts <- counts[, ..sel]
  # Melt
  .m <- melt(counts, id.vars = "FBgn")
  .m[, names(FBgn_list):= lapply(FBgn_list, function(x) as.numeric(FBgn %in% x))]
  .m <- melt(.m, id.vars = c("FBgn", "variable", "value"), variable.name = "cl", value.name = "cl_count")
  # Compute enrichment
  res <- .m[, {
    tab <- table(value>0, cl_count>0)
    if(identical(dim(tab), c(2L,2L)))
    {
      .f <- fisher.test(tab, 
                        alternative = "greater")
      .(OR= .f$estimate,
        pval= .f$p.value)
    }else
      .(OR= as.numeric(NA),
        pval= as.numeric(NA))
  }, .(variable, cl)]
  # padj...
  res[, padj:= p.adjust(pval, method = "fdr"), pval]
  res[, log2OR:= log2(OR)]
  res <- res[, .(variable, cl, log2OR, padj)]
  res[names, variable:= i.name, on= "variable==GO"]
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