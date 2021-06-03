#' GO analysis
#'
#' This function compares GO enrichments for several groups
#'
#' @param FBgn_list A nameds list of FBgn IDs (1 sublist/cluster)
#' @param go_object object containing GO data. see ?vl_fb_go_table_dm6_FB2020_05
#' @param cex cex ballons (usefull to adjust size)
#' @export

vl_GO_clusters <- function(FBgn_list,
                           go_object= vl_fb_go_table_dm6_FB2020_05,
                           cex= 1)
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
    
  # Compute go and total counts
  go_current <- copy(go_object)[,.(GO, name, type= unlist(type), FBgn, Symbol)]
  go_current[, total_FBgn:= length(unique(go_object$FBgn))]
  go_current[, total_go:= length(unique(FBgn)), GO]
  
  # Add cluster counts
  FBgn_DT <- rbindlist(lapply(FBgn_list, 
                              function(x) # Select only cluster uniq FBgns that exist in GO object
                                data.table(FBgn= unique(x[x %in% unique(go_object$FBgn)]))), 
                       idcol = "cluster_names")
  FBgn_DT[, total_cluster:= .N, cluster_names]
  go_current <- go_current[rep(seq(nrow(go_current)), each= length(FBgn_list))] # Expand N clusters
  go_current[, cluster_names:= factor(rep(names(FBgn_list), nrow(go_object)))]
  go_current <- merge(go_current, FBgn_DT, c("cluster_names", "FBgn")) # Keep only lines for which FBgn exist in clusters
  go_current[, go_cluster:= length(unique(FBgn)), .(GO, cluster_names)]
  
  # Compute statistically enriched GOs
  go_current <- go_current[, .(Symbol= .(Symbol), FBgn= .(FBgn)), 
                           .(GO, name, type, cluster_names, total_go, total_FBgn, total_cluster, go_cluster)]
  go_current[, c("estimate", "pvalue"):= {
    fisher.test(matrix(c(go_cluster, # go+ cluster +
                         total_cluster-go_cluster, # go- cluster +
                         total_go-go_cluster, # go+ cluster -
                         total_FBgn-total_cluster-total_go+go_cluster), # go- cluster -
                       nrow=2, 
                       byrow = T))[c("estimate", "p.value")]
  }, total_go:go_cluster]
  go_current[, '-log10(padj)':= -log10(p.adjust(pvalue, method = "fdr"))]
  go_current[, log2OR:= log2(estimate)]
  
  # Generate plot table
  pl <- go_current[`-log10(padj)`>5 & log2OR>0]
  Cc <- colorRamp2(range(pl$`-log10(padj)`), 
                   colors = c("blue", "red"))
  pl[, col:= Cc(`-log10(padj)`)]
  pl[, size:= cex*log2OR]
  pl[, x:= seq(0, 1, length.out = length(unique(pl$cluster_names)))[.GRP], cluster_names]
  setorderv(pl, 
            c("cluster_names", "-log10(padj)", "log2OR", "GO"), 
            order = c(1, -1, -1, 1))
  pl[, y:= as.numeric(min(.I)), GO]
  pl[, y:= seq(1, 0, length.out = length(unique(pl$GO)))[.GRP], keyby= y]
  
  # PLOT init
  par(mai = c(max(strwidth(pl$cluster_names, "inches"))+0.5,
              max(strwidth(pl$name, "inches"))+0.5,
              0.25,
              1),
      xaxs= "i",
      yaxs= "i")
  
  # Lines
  plot.new()
  abline(v= unique(pl$x))
  abline(h= unique(pl$y))
  # Points
  points(pl$x, 
         pl$y, 
         cex= pl$size, 
         col= adjustcolor(pl$col, 0.8),
         pch= 16,
         xpd= T)
  axis(1,
       lty= 0,
       at = unique(pl[, .(cluster_names, x)])$x, 
       labels = unique(pl[, .(cluster_names, x)])$cluster_names, 
       las= 2)
  axis(2,
       lty= 0,
       at = unique(pl[, .(name, y)])$y, 
       labels = unique(pl[, .(name, y)])$name, 
       las= 2)
  # Legend pval
  rasterImage(matrix(Cc(seq(min(pl$`-log10(padj)`), 
                            max(pl$`-log10(padj)`), 
                            length.out = 100)), 
                     ncol= 1), 
              xleft = 1.075,
              ybottom = 0.96,
              xright = 1.125,
              ytop = 0.8, 
              xpd= T)
  text(1.05,
       1-grconvertY(0.5, "line", "ndc"), 
       pos= 4,
       "-log10(padj)", 
       xpd= T, 
       offset= 0, 
       cex= 0.8)
  ticks <- axisTicks(range(pl$`-log10(padj)`), log= F)
  segments(1.125+0.005, 
           0.8+ticks/max(pl$`-log10(padj)`)*(0.96-0.8), 
           1.125+0.01,
           0.8+ticks/max(pl$`-log10(padj)`)*(0.96-0.8),
           xpd= T, 
           lend= 2)
  text(1.125+0.01, 
       0.8+ticks/max(pl$`-log10(padj)`)*(0.96-0.8),
       labels = ticks,
       cex= 0.6,
       offset= 0.1,
       pos= 4,
       xpd= T)
  # Legend balloons
  scale <- axisTicks(range(pl$size), log= F)
  points(rep(1.1, length(scale)),
         seq(0.4, 0.65, length.out = length(scale)), 
         xpd= T,
         col= "black",
         cex= scale,
         pch= 16)
  text(1.05,
       0.71, 
       pos= 4,
       "log2(OR)", 
       xpd= T, 
       offset= 0, 
       cex= 0.8)
  text(1.175,
       seq(0.4, 0.65, length.out = length(scale)), 
       labels = scale,
       xpd= T,
       pos= 4,
       offset= 0)
  
  # RETURN
  invisible(list(data= go_current, 
                 plot= pl, 
                 legend_ticks= ticks, 
                 ballon_ticks= scale))
}