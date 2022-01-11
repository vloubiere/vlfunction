#' Heatmap
#'
#' This function plots a heatmap and returns an S3 object of type vl_heatmap_obj which can be modified and used for later plotting.
#' 
#' @param x Data to plot. Can be one of matrix, data.table, or formula.
#' @param breaks Heatmap breaks, which should be the same length as the color vector
#' @param show_rownames Plot row names? Default= T
#' @param show_colnames Plot col names? Default= T
#' @param col Color scale. default= c("cornflowerblue", "white", "red")
#' @param main Title. Default= NA
#' @param auto_margins Should auto margins be computed and used?
#' @param cluster_rows Should rows be clustered? Default= T
#' @param clustering_distance_rows Method for clustering distance rows. Can be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
#' @param cutree_rows Number of cuts for rows tree. Default= 1L
#' @param show_row_dendrogram Plot row dendrogram?
#' @param show_row_clusters Plot row cluster names/lines?
#' @param cluster_cols Should columns be clustered? Default= T
#' @param clustering_distance_cols Method for clustering distance cols. Can be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
#' @param cutree_cols  Number of cuts for cols tree. default= 1L
#' @param show_col_dendrogram Plot col dendrogram?
#' @param show_col_clusters Plot col cluster names/lines?
#' @param clustering_method Clustering method to use. Must be one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".
#' @param kmeans_k Number of row kmean clusters. If specified, takes over row clustering.
#' @param display_numbers Display numbers on heatmap? Default= F
#' @param display_numbers_FUN Function to apply before displaying numbers on heatmap.
#' @param display_numbers_cex cex display numbers 
#' @param legend_title Character to plot as title legend
#' @param na_col Color for na values. Default= "lightgrey"
#' @param plot Should the heatmap be plotted? Default= T
#' @examples
#' # Create test matrix
#' set.seed(1234)
#' test = matrix(rnorm(200), 20, 10)
#' test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
#' test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
#' test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
#' colnames(test) = paste("Test", 1:10, sep = "")
#' rownames(test) = paste("Gene", 1:20, sep = "")
#' 
#' par(mfrow=c(2,2))
#' res <- vl_heatmap(test, kmeans_k = 4, cutree_rows = 3, cutree_cols = 5, legend_title = "log2FoldChange", main= "test")
#' res$result_DT[]
#' vl_heatmap(test, cutree_rows = 3, cutree_cols = 2, show_col_clusters = F, legend_title = "log2FoldChange", main= "test")
#' vl_heatmap(test, cutree_rows = 3, cutree_cols = 2, show_col_dendrogram = F, legend_title = "log2FoldChange", main= "test")
#' obj <- vl_heatmap(test, cutree_rows = 3, cutree_cols = 2, show_row_clusters = F, show_row_dendrogram = F,show_col_clusters = F, show_col_dendrogram = F, legend_title = "log2FoldChange", main= "test", plot= F)
#' plot(obj, cutree_rows = 3, cutree_cols = 2, show_row_clusters = F, show_row_dendrogram = F,show_col_clusters = F, show_col_dendrogram = F, legend_title = "log2FoldChange", main= "test")
#' 
#' @return An object of class vl_heatmap_obj containing clustering data and can be used for plotting using plot(obj, ...)
#' @export
vl_heatmap <- function(x, ...) UseMethod("vl_heatmap")

#' @describeIn vl_heatmap Just a wrapper that collapses DT to a suitable matrix
#' @export
vl_heatmap.data.table <- function(x, rownames= names(x)[1], ...)
{
  if(rownames %in% names(x))
  {
    setcolorder(x, rownames)
    mat <- as.matrix(x, 1)
  }else
    mat <- as.matrix(x)
  vl_heatmap.matrix(mat, ...)
}

#' @describeIn vl_heatmap Core function that clusters input matrix and plots a heatmap
#' @export
vl_heatmap.matrix <- function(x,
                              breaks= seq(min(x, na.rm= T), max(x, na.rm= T), length.out= length(col)),
                              show_rownames= T,
                              show_colnames= T,
                              col= c("cornflowerblue", "white", "red"),
                              main= NA,
                              auto_margins= T,
                              cluster_rows= T,
                              clustering_distance_rows= "euclidean",
                              cutree_rows = 1,
                              show_row_dendrogram= T,
                              show_row_clusters= T,
                              cluster_cols = T,
                              clustering_distance_cols = "euclidean",
                              cutree_cols = 1,
                              show_col_dendrogram= T,
                              show_col_clusters= T,
                              clustering_method = "complete",
                              kmeans_k= NA,
                              display_numbers= F,
                              display_numbers_FUN= function(x) round(x, 2),
                              display_numbers_cex= 0.5,
                              legend_title= NA,
                              na_col= "lightgrey",
                              plot= T)
{
  if(is.null(rownames(x)))
    rownames(x) <- seq(nrow(x))
  #------------------------####
  # Init informative result DT
  #------------------------####
  clustered_x <- x
  result_DT <- as.data.table(x, keep.rownames = T)
  result_DT <- melt(result_DT, id.vars = "rn")
  
  #------------------------####
  # Clustering rows
  #------------------------####
  if(cluster_rows & nrow(clustered_x)>1)
  {
    set.seed(3453)
    if(missing(kmeans_k))
    {
      if(clustering_distance_rows %in% c("pearson", "spearman"))
        .d <- as.dist(1 - cor(t(clustered_x), 
                              use= "pairwise.complete.obs", 
                              method= clustering_distance_rows)) else
                                .d <- dist(clustered_x, method = clustering_distance_rows)
                              # Hierarchical clustering
                              rcl <- hclust(.d, method = clustering_method) 
                              clustered_x <- clustered_x[rcl$order,]
                              # Add to DT
                              result_DT[data.table(rcl$labels, 
                                                   cutree(rcl, k = cutree_rows)), rn_cl:= i.V2, on= "rn==V1"]
                              # Extract dend
                              rdend <- data.table::as.data.table(ggdendro::dendro_data(rcl, 
                                                                                       type = "rectangle", 
                                                                                       rotate= T)$segments)
                              
    }else
    {
      # Kmeans clustering
      rcl <- kmeans(clustered_x, centers = kmeans_k) 
      clustered_x <- clustered_x[order(rcl$cluster),]
      # Add to DT
      result_DT[data.table(names(rcl$cluster), 
                           rcl$cluster), rn_cl:= i.V2, on= "rn==V1"]
    }
  }
  
  #------------------------####
  # Clustering cols
  #------------------------####
  if(cluster_cols & ncol(clustered_x)>1)
  {
    set.seed(3453)
    if(clustering_distance_cols %in% c("pearson", "spearman"))
      .d <- as.dist(1 - cor(clustered_x, 
                            use= "pairwise.complete.obs", 
                            method= clustering_distance_cols)) else
                              .d <- dist(t(clustered_x), method = clustering_distance_cols)
                            # Hierarchical clustering
                            ccl <- hclust(.d, 
                                          method = clustering_method)
                            clustered_x <- clustered_x[,ccl$order]
                            # Extract dend
                            cdend <- data.table::as.data.table(ggdendro::dendro_data(ccl,
                                                                                     type = "rectangle", 
                                                                                     rotate= T)$segments)
                            # Add to DT
                            result_DT[data.table(ccl$labels, 
                                                 cutree(ccl, k = cutree_cols)), var_cl:= i.V2, on= "variable==V1"]
  }
  
  #------------------------####
  # PLOT
  #------------------------####
  obj <- as.list(environment())
  if(plot)
  {
    class(obj) <- "vl_heatmap_pl"
    plot(obj)
  }
  
  class(obj) <- "vl_heatmap_obj"
  invisible(obj)
}

#' @export
plot.vl_heatmap_obj <- function(obj,
                                add,
                                add_inherit_row_order= T,
                                add_inherit_col_order= F,
                                breaks= seq(min(obj$x, na.rm= T), max(obj$x, na.rm= T), length.out= length(col)),
                                show_rownames= T,
                                show_colnames= T,
                                col= c("cornflowerblue", "white", "red"),
                                main= NA,
                                auto_margins= T,
                                cutree_rows = 1,
                                show_row_dendrogram= T,
                                show_row_clusters= T,
                                cutree_cols = 1,
                                show_col_dendrogram= T,
                                show_col_clusters= T,
                                display_numbers= F,
                                display_numbers_FUN= function(x) round(x, 2),
                                display_numbers_cex= 0.5,
                                legend_title= NA)
{
  if(!missing(add))
  {
    if(!is.matrix(add))
      stop("add should be a matrix with either same number of rows and/or cols as original input") else
        obj$x <- add
    if(add_inherit_row_order)
      if(nrow(add)!=nrow(obj$x))
        stop("add_inherit_row_order is TRUE but nrow(add)!=nrow(x)") else if(obj$cluster_rows)
        {
          if(class(obj$rcl)=="hclust")
            add <- add[obj$rcl$order,, drop=F] else if(class(obj$rcl)=="kmeans")
              add <- add[order(obj$rcl$cluster),, drop=F]
        }
    if(add_inherit_col_order)
      if(ncol(add)!=ncol(obj$x))
        stop("add_inherit_col_order is TRUE but ncol(add)!=ncol(x)") else if(obj$cluster_cols)
          add <- add[,obj$ccl$order, drop=F]
        show_row_dendrogram <- F
        show_col_dendrogram <- F
        obj$clustered_x <- add
  }
  args <- ls()
  args <- args[args %in% names(obj)]
  for(arg in args)
    obj[[arg]] <- get(arg)
  class(obj) <- "vl_heatmap_pl"
  plot(obj)
}

# Default plotting function
plot.vl_heatmap_pl <- function(obj)
{
  list2env(obj, envir = environment())
  # Margins
  if(auto_margins)
  {
    bot <- 1
    left <- 1
    top <- 1.5
    right <- 0.5
    if(show_colnames)
      bot <- bot+grconvertY(max(strwidth(colnames(clustered_x), "inches")), "inches", "lines")
    if(show_rownames)
      left <- left+grconvertX(max(strwidth(rownames(clustered_x), "inches")), "inches", "lines")
    if(is.character(main))
      top <- top+grconvertY(strheight(main, 
                                      cex = 0.8, 
                                      units = "inches"), "inches", "lines")
    if(exists("cdend"))
      if(show_col_dendrogram)
        top <- top+2 else if(show_col_clusters)
          top <- top+1
    if(exists("rdend") & show_row_dendrogram)
      right <- right+2
    if(show_row_clusters & (is.numeric(kmeans_k) | exists("rdend")))
      right <- right+1
    leg_width <- grconvertX(strwidth(legend_title, units = "inches"), "inches", "lines")
    if(leg_width>3.5)
      right <- right+leg_width else
        right <- right+2.5
    par(mar= c(bot, left, top, right))
  }
  
  # Image
  Cc <- circlize::colorRamp2(breaks, 
                             colors= col)
  im <- clustered_x
  im[!is.na(im)] <- Cc(clustered_x[!is.na(clustered_x)])
  im[is.na(im)] <- na_col
  plot.new()
  plot.window(xlim = c(0.5,ncol(im)+0.5),
              ylim = c(0.5,nrow(im)+0.5),
              xaxs= "i",
              yaxs= "i")
  rasterImage(im,
              xleft = 0.5,
              ybottom = 0.5,
              xright = ncol(im)+0.5,
              ytop = nrow(im)+0.5,
              interpolate = F)
  
  # Plot numbers
  if(display_numbers)
    text(c(col(clustered_x)),
         rev(c(row(clustered_x))),
         display_numbers_FUN(c(clustered_x)),
         cex= display_numbers_cex,
         offset= 0)
  
  # Plot dendro and cuts
  adjpos <- function(x, ax)
  {
    .f <- switch(ax, 
                 "x"= graphics::grconvertX,
                 "y"= graphics::grconvertY)
    x <- x/max(x)*.f(2, "lines", "ndc") # Norm
    x <- x+.f(switch(ax, 
                     "x"= par("usr")[2],
                     "y"= par("usr")[4]), "user", "ndc") # Position
    x <- .f(x, "ndc", "user") # To user
    return(x)
  }
  if(exists("rdend"))
  {
    if(show_row_dendrogram)
      segments(adjpos(rdend$y, ax= "x"),
               -rdend$x+nrow(im)+1,
               adjpos(rdend$yend, ax= "x"),
               -rdend$xend+nrow(im)+1,
               xpd= T)
    cuts <- data.table(rev(cutree(rcl, cutree_rows)[rcl$order]))
    cuts[, V2:= .I]
    last_cut <- cuts[.N, V1]
    cuts[, {
      if(V1!=last_cut)
        abline(h= max(V2+0.5))
      if(cutree_rows>1 & show_row_clusters)
        text(grconvertX(1, "npc", "user"),
             mean(V2), 
             V1[1],
             pos= 4,
             xpd= T,
             cex= 1.5, 
             offset = ifelse(show_row_dendrogram, 1, 0.25))
    }, V1]
  }
  if(exists("cdend"))
  {
    if(show_col_dendrogram)
      segments(cdend$x,
               adjpos(cdend$y, "y"),
               cdend$xend,
               adjpos(cdend$yend, "y"),
               xpd= T)
    cuts <- data.table(cutree(ccl, cutree_cols)[ccl$order])
    cuts[, V2:= .I]
    last_cut <- cuts[.N, V1]
    cuts[, {
      if(V1!=last_cut)
        abline(v= max(V2+0.5))
      if(cutree_cols>1 & show_col_clusters)
        text(mean(V2),
             grconvertY(1, "npc", "user"),
             V1[1],
             pos= 3,
             xpd= T,
             cex= 1.5, 
             offset = ifelse(show_col_dendrogram, 1, 0.25))
    }, V1]
  }
  if(is.numeric(kmeans_k) & show_row_clusters)
  {
    cuts <- data.table(rcl$cluster)
    cuts[order(V1, decreasing = T), V2:= .I]
    last_cut <- cuts[.N, V1]
    cuts[, {
      if(V1>1)
        abline(h= max(V2+0.5))
      text(grconvertX(1, "npc", "user"),
           mean(V2),
           V1[1],
           pos= 4,
           xpd= T,
           cex= 1.5,
           offset = 0.25)
    }, V1]
  }
  
  # Plot axes
  if(show_colnames & !is.null(colnames(clustered_x)))
    axis(1,
         at= seq(ncol(clustered_x)),
         labels = colnames(clustered_x),
         lwd= 0,
         las= 2,
         line= -0.5)
  if(show_rownames & !is.null(rownames(clustered_x)))
    axis(2,
         at= seq(nrow(clustered_x)),
         labels = rev(rownames(clustered_x)),
         lwd= 0,
         las= 2,
         line= -0.5)
  
  # Plot legend
  xleft <- grconvertX(1, "npc", "ndc")+grconvertX(0.5, "lines", "ndc")
  if(exists("rdend") & show_row_dendrogram)
    xleft <- xleft+grconvertX(2, "lines", "ndc")
  if(show_row_clusters & (is.numeric(kmeans_k) | exists("rdend")))
    xleft <- xleft+grconvertX(1, "lines", "ndc")
  xright <- xleft+grconvertX(1, "lines", "ndc")
  xleft <- grconvertX(xleft, "ndc", "user")
  xright <- grconvertX(xright, "ndc", "user")
  ytop <- grconvertY(1, "npc", "ndc")-grconvertY(1.25, "lines", "ndc")
  yleg <- ytop+grconvertY(1, "lines", "ndc")
  yleg <- grconvertY(yleg, "ndc", "user")
  ybottom <- ytop-grconvertY(5, "lines", "ndc")
  ytop <- grconvertY(ytop, "ndc", "user")
  ybottom <- grconvertY(ybottom, "ndc", "user")
  rasterImage(matrix(rev(Cc(seq(min(breaks), max(breaks), length.out = 101)))),
              xleft,
              ybottom,
              xright,
              ytop,
              xpd=T)
  ticks <- axisTicks(range(breaks), log=F)
  ymin.ticks <- ybottom+(min(ticks)-min(breaks))/diff(range(breaks))*(ytop-ybottom)
  ymax.ticks <- ybottom+(max(ticks)-min(breaks))/diff(range(breaks))*(ytop-ybottom)
  text(xright,
       seq(ymin.ticks, ymax.ticks, length.out = length(ticks)),
       labels = ticks,
       pos=4,
       xpd= T,
       cex= 0.6,
       offset= 0.25)
  text(xleft,
       yleg,
       labels = legend_title,
       pos=4,
       xpd= T,
       cex= 0.8,
       offset= 0)
  text(grconvertX(0.5, "npc", "user"),
       grconvertY(grconvertY(1, "nfc", "ndc")-grconvertY(1, "lines", "ndc"), "ndc", "user"),
       cex= 1.25,
       main, 
       xpd= T)
  box(lwd= 0.25)
}
