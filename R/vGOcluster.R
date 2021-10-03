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
                           cex= 1,
                           padj_cutoff= 1e-5,
                           N_top= Inf,
                           auto_margin= T)
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
  Cc <- colorRamp2(pval_lims, 
                   colors = c("blue", "red"))
  pl[, col:= Cc(`-log10(pval)`)]
  # Points scaling factor
  pl[, size:= (2^cor_log2OR)]
  pl[, size:= cex*log2(size+1)]
  size_lims <- range(pl$size)
  # Y coordinates
  setorderv(pl, 
            c("cluster_name", "-log10(pval)", "cor_log2OR", "GO"), 
            order = c(1, -1, -1, 1))
  pl[, y:= as.numeric(min(.I)), GO]
  pl[, y:= seq(1, 0, length.out = length(unique(pl$GO)))[.GRP], keyby= y]
  # X coordinates
  pl[, x:= seq(0, 1, length.out = length(unique(pl$cluster_name)))[.GRP], keyby= cluster_name]
  
  #-----------------------------#
  # PLOT
  #-----------------------------#
  if(auto_margin)
    par(mai = c(max(strwidth(pl$cluster_name, "inches"))+0.5,
                max(strwidth(pl$name, "inches"))+0.5,
                0.5,
                strwidth("-log10(pval)", "inches")+0.25),
        xaxs= "i",
        yaxs= "i")

  # Lines
  plot.new()
  segments(unique(pl$x),
           par("usr")[3],
           unique(pl$x),
           par("usr")[4])
  segments(par("usr")[1],
           unique(pl$y),
           par("usr")[2],
           unique(pl$y))
  # Points
  points(pl$x, 
         pl$y, 
         cex= pl$size, 
         col= adjustcolor(pl$col, 0.8),
         pch= 16,
         xpd= T)
  axis(1,
       lty= 0,
       at = unique(pl$x), 
       labels = unique(pl$cluster_name), 
       las= 2)
  axis(2,
       lty= 0,
       at = unique(pl$y), 
       labels = unique(pl$name),
       las= 2)
  # Legend pval
  xleft <- grconvertX(1, "npc", "inches")+grconvertX(1, "lines", "inches")
  xright <- xleft+grconvertX(1, "lines", "inches")
  xleft <- grconvertX(xleft, "inches", "npc")
  xright <- grconvertX(xright, "inches", "npc")
  ybottom <- 0.7
  ytop <- grconvertY(1, "npc", "inches")-grconvertY(1, "chars", "inches")
  ytop <- grconvertY(ytop, "inches", "npc")
  rasterImage(matrix(rev(Cc(seq(min(pval_lims), max(pval_lims), length.out = 101)))),
              xleft,
              ybottom,
              xright,
              ytop,
              xpd=T)
  ticks <- axisTicks(pval_lims, log=F)
  ymin.ticks <- ybottom+(min(ticks)-pval_lims[1])/diff(pval_lims)*(ytop-ybottom)
  ymax.ticks <- ybottom+(max(ticks)-pval_lims[1])/diff(pval_lims)*(ytop-ybottom)
  text(xright,
       seq(ymin.ticks, ymax.ticks, length.out = length(ticks)),
       labels = ticks,
       pos=4,
       xpd= T,
       cex= 0.6,
       offset= 0.25)
  text(xleft,
       1-0.5*strheight("A"),
       labels = "-log10(pval)",
       pos=4,
       xpd= T,
       cex= 0.8,
       offset= 0)
  
  # Legend balloons
  scale <- axisTicks(size_lims, log= F)
  maxBalloonInch <- strheight("A", units = "inches", cex= max(scale))*0.75
  bx <- grconvertX(xleft, "npc", "inches")+maxBalloonInch/2
  bx <- grconvertX(bx, "inches", "npc")
  btop <- grconvertY(0.6, "npc", "inches")-maxBalloonInch/2
  bbot <- btop-maxBalloonInch*(length(scale)-1)
  btop <- grconvertY(btop, "inches", "npc")
  bbot <- grconvertY(bbot, "inches", "npc")
  points(rep(bx, length(scale)),
         seq(bbot, btop, length.out = length(scale)), 
         xpd= T,
         col= "black",
         cex= scale,
         pch= 16)
  bx <- grconvertX(bx, "npc", "inches")+maxBalloonInch/2
  bx <- grconvertX(bx, "inches", "npc")
  text(rep(bx, length(scale)),
       seq(bbot, btop, length.out = length(scale)),
       labels= scale,
       pos= 4,
       xpd= T,
       offset= 0,
       cex= 0.8)
  tity <- grconvertY(0.6, "npc", "inches")+grconvertY(1, "chars", "inches")/2
  tity <- grconvertY(tity, "inches", "npc")
  text(xleft,
       tity,
       labels = "log2OR",
       pos= 4,
       xpd= T,
       cex= 0.8,
       offset= 0)
  
  # RETURN
  invisible(list(data= res, 
                 plot= pl))
}
