#' GO analysis 
#'
#' This function compute GO enrichment in a group of genes
#'
#' @param FBgn_list A vector of FBgn IDs
#' @param go_object object containing GO data. see ?vl_fb_go_table_dm6_FB2020_05
#' @param FBgn_universe Vector of FBgn IDs to which the analysis will be restricted. default= "all" (i.e. all the FBGns present in the table)
#' @param main title
#' @export

vl_GO_enrich <- function(FBgn_vector,
                         go_object= vl_fb_go_table_dm6_FB2020_05,
                         FBgn_universe= "all",
                         go_type= "all",
                         main= "")
{
  # Checks
  if(!is.vector(FBgn_vector))
    stop(paste0("FBgn_vector should be a vector FBgn IDs"))
  if(!all(grepl("FBgn", FBgn_vector)))
    stop(paste0("Some FBgn_vector components are not FBgn IDs"))
  if(length(FBgn_universe)==1)
  {
    if(FBgn_universe=="all")
      FBgn_universe <- unique(go_object$FBgn)
  }
  if(length(go_type)!=1)
    stop("length(go_type)!=1")
  if(!(go_type %in% c("all", "biological_process", "cellular_component", "molecular_function")))
     stop("go_type should be one of 'all', 'biological_process', 'cellular_component', 'molecular_function'")
    
  # Restrict GO object to universe
  go_current <- copy(go_object)[FBgn %in% FBgn_universe, .(GO, name, type= unlist(type), FBgn, Symbol)]
  
  # Restrict to type
  if(go_type!="all")
    go_current <- go_current[type==go_type]
  
  # Compute go and total counts
  go_current[, total_FBgn:= length(FBgn_universe)]
  go_current[, total_go:= length(unique(FBgn)), GO]
  
  # Add cluster counts
  go_current <- go_current[FBgn %in% FBgn_vector]
  go_current[, total_cluster:= length(unique(FBgn))]
  go_current[, go_cluster:= length(unique(FBgn)), GO]
  
  # Compute statistically enriched GOs
  go_current <- go_current[, .(Symbol= .(Symbol), FBgn= .(FBgn)), 
                           .(GO, name, type, total_go, total_FBgn, total_cluster, go_cluster)]
  go_current[, c("estimate", "pvalue"):= {
    fisher.test(matrix(c(go_cluster, # go+ cluster +
                         total_cluster-go_cluster, # go- cluster +
                         total_go-go_cluster, # go+ cluster -
                         total_FBgn-total_cluster-total_go+go_cluster), # go- cluster -
                       nrow=2, 
                       byrow = T), 
                alternative = "greater")[c("estimate", "p.value")]
  }, total_go:go_cluster]
  go_current[, '-log10(padj)':= -log10(p.adjust(pvalue, method = "fdr"))]
  go_current[, log2OR:= log2(estimate)]
  
  #----------------------------#
  # Generate plot table
  #----------------------------#
  pl <- go_current[`-log10(padj)`>5 & log2OR>0]
  # Handle infinite values
  pl[, cor_log2OR:= log2OR]
  pl[log2OR==Inf, cor_log2OR:= max(pl$log2OR[is.finite(pl$log2OR)], na.rm= T)]
  Cc <- colorRamp2(range(pl$cor_log2OR, na.rm = T), 
                   colors = c("blue", "red"))
  pl[, col:= Cc(log2OR)]
  setorderv(pl, "-log10(padj)")
  
  #----------------------------#
  # PLOT
  #----------------------------#
  par(mai= c(1, 
             max(strwidth(pl$name, "inches"))+1, 
             0.5, 
             1))
  barplot(pl$`-log10(padj)`, 
          col= pl$col, 
          horiz = T,
          xlab= "-log10(padj)", 
          names.arg = pl$name,
          las= 1, 
          main= main, 
          space = 0,
          border= "white")
  rasterImage(matrix(Cc(seq(min(pl$cor_log2OR), 
                            max(pl$cor_log2OR), 
                            length.out = 100)), 
                     ncol= 1), 
              xleft = grconvertX(1.075, "npc", "user"),
              ybottom = grconvertY(0.93, "npc", "user"),
              xright = grconvertX(1.125, "npc", "user"),
              ytop = grconvertY(0.75, "npc", "user"), 
              xpd= T)
  text(grconvertX(1.075, "npc", "user"),
       grconvertY(0.955, "npc", "user"), 
       pos= 4,
       "log2(OR)", 
       xpd= T, 
       offset= 0, 
       cex= 0.8)
  ticks <- axisTicks(range(pl$cor_log2OR), log= F)
  at <- grconvertY(0.75+(ticks-min(pl$cor_log2OR))/(max(pl$cor_log2OR)-min(pl$cor_log2OR))*(0.93-0.75), "npc", "user")
  segments(grconvertX(1.125+0.005, "npc", "user"),
           at,
           grconvertX(1.125+0.01, "npc", "user"),
           at,
           xpd= T,
           lend= 2)

  text(grconvertX(1.125+0.01, "npc", "user"), 
       at,
       labels = ticks,
       cex= 0.6,
       offset= 0.1,
       pos= 4,
       xpd= T)
  
  # RETURN
  invisible(list(data= go_current, 
                 plot= pl, 
                 legend_ticks= ticks,
                 ticks_at= at))  
  
}

#' GO analysis
#'
#' This function compares GO enrichments for several groups
#'
#' @param FBgn_list A nameds list of FBgn IDs (1 sublist/cluster)
#' @param go_object object containing GO data. see ?vl_fb_go_table_dm6_FB2020_05
#' @param all_FBgns Vector of FBgns to be used in the universe, typically all the genes tested with DESeq. Default= "all" FBgns present in go_object.
#' @param go_type The type of go to be considered. Default "all" means "biological_process" AND "cellular_component" AND "molecular_function". 
#' @param cex cex ballons (usefull to adjust size)
#' @param padj_cutoff padjust cutoff applied to GOs
#' @param N_top Select only N top GOs/cluster
#' @param auto_margin Compute and apply optimal margins. Default= T
#' @export

vl_GO_clusters <- function(FBgn_list,
                           go_object= vl_fb_go_table_dm6_FB2020_05,
                           all_FBgns= "all",
                           go_type= "all",
                           padj_cutoff= 1e-5,
                           N_top= Inf,
                           auto_margins= T,
                           cex.balloons= 4,
                           main= NA)
{
  # Checks
  if(!is.list(FBgn_list))
    stop(paste0("FBgn_list should be a named list of gene FBgn IDs"))
  if(is.null(names(FBgn_list)))
    names(FBgn_list) <- paste0("Group_", seq(FBgn_list))
  if(length(unique(names(FBgn_list)))!=length(FBgn_list))
    stop("All FBgn list names should be unique!")
  if(!all(grepl("FBgn", unlist(FBgn_list))))
    stop(paste0("Some FBgn_list components are not FBgn IDs"))
  if(length(all_FBgns)==1)
  {
    if(all_FBgns!="all")
      stop(paste0("When length(all_FBgns)==1, only valid value is 'all'"))
    else
      all_FBgns <- unique(go_object$FBgn)
  }
  if(length(go_type)!=1 |
     !(go_type %in% c("all", "biological_process", "cellular_component", "molecular_function")))
    stop("go_type should be one of 'all', 'biological_process', 'cellular_component', 'molecular_function'")
  
  # Restrict GO object to universe
  go_current <- copy(go_object)[FBgn %in% all_FBgns, .(GO, name, type= unlist(type), FBgn, Symbol)]
  perc_noGO <- length(which(!unlist(FBgn_list) %in% go_current$FBgn))/length(unlist(FBgn_list))*100
  perc_noGO <- round(perc_noGO, 1)
  print(paste0(perc_noGO, "% of provided FBgns were not found in GO database!\n"))
  
  # Restrict to type
  if(go_type!="all")
    go_current <- go_current[type==go_type]
  
  # Compute go and total counts
  tab <- copy(go_current)
  tab[, total_FBgn:= length(unique(go_current$FBgn))]
  tab[, total_go:= length(unique(FBgn)), GO]
  
  # Add cluster counts
  FBgn_DT <- rbindlist(lapply(FBgn_list, as.data.table), idcol = "cluster_name")
  FBgn_DT[, total_cluster:= .N, cluster_name]
  tab <- merge(tab,
               FBgn_DT, 
               by.x= "FBgn",
               by.y= "V1")
  tab[, go_cluster:= length(unique(FBgn)), .(GO, cluster_name)]
  
  # Compute statistically enriched GOs
  res <- tab[, .(Symbol= .(Symbol), FBgn= .(FBgn)), 
             .(GO, name, type, cluster_name, total_go, total_FBgn, total_cluster, go_cluster)]
  res[, c("estimate", "pvalue"):= {
    fisher.test(matrix(c(go_cluster, # go+ cluster +
                         total_cluster-go_cluster, # go- cluster +
                         total_go-go_cluster, # go+ cluster -
                         total_FBgn-total_cluster-total_go+go_cluster), # go- cluster -
                       nrow=2, 
                       byrow = T), 
                alternative = "greater")[c("estimate", "p.value")]
  }, total_go:go_cluster]
  res[, "-log10(pval)":= -log10(pvalue)]
  res[, padj:= p.adjust(pvalue, method = "fdr")]
  res[, log2OR:= log2(estimate)]
  
  #----------------------------------#
  # Generate plot table
  #----------------------------------#
  sel <- res[any(padj<=padj_cutoff & log2OR>0), GO]
  pl <- res[GO %in% sel & padj<0.05 & log2OR>0]
  # Select top GO/cluster
  setorderv(pl, "-log10(pval)", order = -1)
  pl <- pl[GO %in% pl[, GO[seq(.N)<=N_top], cluster_name]$V1]
  # Make cluster names as factor
  pl[, cluster_name:= factor(cluster_name, levels = names(FBgn_list))]
  # Handle infinite values
  pl[, cor_log2OR:= log2OR]
  pl[log2OR==Inf, cor_log2OR:= max(pl$log2OR[is.finite(pl$log2OR)], na.rm= T)]
  # Colors
  pval_lims <- range(pl$`-log10(pval)`)
  # Y ordering
  setorderv(pl, 
            c("cluster_name", "-log10(pval)", "cor_log2OR", "GO"), 
            order = c(1, -1, -1, 1))
  pl[, name:= factor(name, levels= unique(pl$name))]
  
  #-----------------------------#
  # PLOT
  #-----------------------------#
  x <- as.matrix(dcast(pl, name~cluster_name, value.var = "cor_log2OR"), 1)
  color_var <- as.matrix(dcast(pl, name~cluster_name, value.var = "-log10(pval)"), 1)
  vl_balloons_plot(x = x,
                   color_var= color_var,
                   col= c("blue", "red"),
                   main= main,
                   cex.balloons = cex.balloons,
                   auto_margins = auto_margins,
                   balloon_size_legend= "OR (log2)",
                   balloon_col_legend = "padj (-log10)")
}

