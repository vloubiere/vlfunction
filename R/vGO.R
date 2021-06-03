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
                         main= "")
{
  # Checks
  if(!is.vector(FBgn_vector))
    stop(paste0("FBgn_vector should be a vector FBgn IDs"))
  if(!all(grepl("FBgn", FBgn_vector)))
    stop(paste0("Some FBgn_vector components are not FBgn IDs"))
  if(length(FBgn_universe)==1)
    if(FBgn_universe=="all")
      FBgn_universe <- unique(go_object$FBgn)  
    
    # Restrict GO object to universe
    go_current <- copy(go_object)[FBgn %in% FBgn_universe, .(GO, name, type= unlist(type), FBgn, Symbol)]
    
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
                         byrow = T))[c("estimate", "p.value")]
    }, total_go:go_cluster]
    go_current[, '-log10(padj)':= -log10(p.adjust(pvalue, method = "fdr"))]
    go_current[, log2OR:= log2(estimate)]
    
    # Generate plot table
    pl <- go_current[`-log10(padj)`>5 & log2OR>0]
    Cc <- colorRamp2(range(pl$log2OR), 
                     colors = c("blue", "red"))
    pl[, col:= Cc(log2OR)]
    setorderv(pl, "-log10(padj)")
    
    # PLOT
    par(mai= c(1, 
               max(strwidth(pl$name, "inches")), 
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
    rasterImage(matrix(Cc(seq(min(pl$log2OR), 
                              max(pl$log2OR), 
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
    ticks <- axisTicks(range(pl$log2OR), log= F)
    segments(grconvertX(1.125+0.005, "npc", "user"), 
             grconvertY(0.75+ticks/max(pl$log2OR)*(0.93-0.75), "npc", "user"), 
             grconvertX(1.125+0.01, "npc", "user"),
             grconvertY(0.75+ticks/max(pl$log2OR)*(0.93-0.75), "npc", "user"),
             xpd= T, 
             lend= 2)
    text(grconvertX(1.125+0.01, "npc", "user"), 
         grconvertY(0.75+ticks/max(pl$log2OR)*(0.93-0.75), "npc", "user"),
         labels = ticks,
         cex= 0.6,
         offset= 0.1,
         pos= 4,
         xpd= T)
    
    # RETURN
    invisible(list(data= go_current, 
                   plot= pl, 
                   legend_ticks= ticks))
}