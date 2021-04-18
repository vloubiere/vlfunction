#' Heatmap
#'
#' This function plots a violin plot from either matrix,
#' data.table, data.frame or using formula.
#'
#' @param mat Matrix to plot. Row names and col names will be reported in return
#' @param newdata If specified, clusters using matrix and plots new data instead
#' @param breaks Heatmap breaks, which should be the same length as the color vector
#' @param cluster_rows Should matrix rows be clustered? default is T
#' @param clustering_distance_rows Method for clustering distance rows. Can be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
#' @param cutree_rows Number of cuts for rows tree. default= 1z
#' @param cluster_cols Should matrix cols be clustered? default is T
#' @param clustering_distance_cols Method for clustering distance cols. Can be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
#' @param cutree_cols Number of cuts for cols tree. default= 1z
#' @param clustering_method Clustering method to use. Must be one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".
#' @param kmeans_k Number of row kmean clusters. If specified, takes over row clustering.
#' @param display_numbers Display numbers on heatmap? Default= F
#' @param display_numbers_FUN Function to apply before displaying numbers on heatmap.
#' @param display_numbers_cex cex display numbers
#' @param legend_title Character to plot as title legend
#' @param legend_pos vector(bottom, left, top, right). default= c(0.8, 1.15, 1, 1.2),
#' @param scale Scale matrix before clustering? Can be 'none' (default), 'column' or 'row'.
#' @return A data.table containing clustering data!
#' @export

# # Example matrix
# set.seed(123)
# nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
# nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3
# mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
#                   matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
#                   matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
#             rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
#                   matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
#                   matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
#             rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
#                   matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
#                   matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3)))
# mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
# rownames(mat) = paste0("row", seq_len(nr))
# colnames(mat) = paste0("column", seq_len(nc))

vl_heatmap <- function(mat,
                       newdata= NULL,
                       breaks= NULL,
                       show_rownames= T,
                       show_colnames= T,
                       col= c("cornflowerblue", "white", "red"),
                       cluster_rows= T,
                       clustering_distance_rows= "euclidean",
                       cutree_rows = 1,
                       cluster_cols = T,
                       clustering_distance_cols = "euclidean",
                       cutree_cols = 1,
                       clustering_method = "complete",
                       kmeans_k = NA,
                       display_numbers= F,
                       display_numbers_FUN= function(x) round(x, 2),
                       display_numbers_cex= 0.5,
                       legend_title= "legend",
                       legend_pos= c(0.8, 1.15, 1, 1.2),
                       scale = "none")
{
  if(!is.matrix(mat))
    stop("mat must be a matrix object")
  if(is.null(colnames(mat)))
    colnames(mat) <- paste0("column", seq(ncol(mat)))
  if(is.null(rownames(mat)))
    rownames(mat) <- paste0("row", seq(nrow(mat)))
  if(length(unique(rownames(mat)))!=nrow(mat))
    rownames(mat) <- paste0(rownames(mat), "__", seq(nrow(mat)))
  if(!is.character(class(colnames(mat))))
    stop("matrix colnames should be of class character")
  if(!is.character(class(rownames(mat))))
    stop("matrix rownames should be of class character")
  if(!is.null(newdata))
  {
    if(!is.matrix(newdata) | nrow(mat) != nrow(newdata))
      stop("newdata must be a matrix object width same number of rows as mat")
    if(scale!="none")
      warning("Scaling only applied to the clustering data, not newdata!")
  }
  if(!is.null(breaks) & length(breaks)!=length(col))
  {
    col <- colorRampPalette(col)(length(breaks))
    warning("color vector has been interpolated to match breaks vector length")
  }

  #------------------------####
  # Format object
  #------------------------####
  DT <- data.table::as.data.table(mat, keep.rownames = "row")
  DT <- data.table::melt.data.table(DT, measure.vars = colnames(DT)[-1], variable.name = "col")
  if(scale=="row")
    DT[, value:= scale(value), row]
  if(scale=="column")
    DT[, value:= scale(value), col]

  #------------------------####
  # Clustering rows
  #------------------------####
  if(cluster_rows)
  {
    set.seed(3453)
    if(!is.na(kmeans_k))
    {
      rcl <- kmeans(mat, centers = kmeans_k) # Kmeans clustering
      rcl <- data.table::data.table(names(rcl$cluster), rcl$cluster)
      DT[rcl, rcl:= i.V2, on= "row==V1"] # Cluster= kmeans cluster
      DT[, y:= .GRP, keyby= .(rcl, row)] # y coor defined by 1/ cluster 2/rownames
    }else
    {
      if(clustering_distance_rows %in% c("pearson", "spearman"))
      {
        .d <- as.dist(1 - cor(t(mat), use= "pairwise.complete.obs", method= clustering_distance_rows))
      }else
      {
        .d <- dist(mat, method = clustering_distance_rows)
      }
      rcl <- hclust(.d, method = clustering_method) # Hierarchical clustering
      rdend <- data.table::as.data.table(ggdendro::dendro_data(rcl, type = "rectangle", rotate= T)$segments) # Extract dendrogram
      ry <- data.table::data.table(rcl$labels, order(rcl$order))
      DT[ry, c("rcl", "y"):= .(1, i.V2), on= "row==V1"]
      if(!is.na(cutree_rows))
      {
        rcl <- data.table::data.table(names(cutree(rcl, cutree_rows)), cutree(rcl, cutree_rows))
        DT[rcl, rcl:= i.V2, on= "row==V1"]
      }
    }
  }else
  {
    DT[, c("rcl", "y"):= .(1, .GRP), row] # No clustering
  }
  DT[, y:= max(y)-y+1] # Invert value to make them correct for image(x,y,z)
  
  #------------------------####
  # Replace DT with newdata following row clustering
  #------------------------####
  if(!is.null(newdata))
  {
    # Rename dims
    rownames(newdata) <- rownames(mat)
    if(is.null(colnames(newdata)))
      colnames(newdata) <- paste0("column", seq(ncol(newdata)))
    # Format
    add <- data.table::as.data.table(newdata, 
                                     keep.rownames = "row")
    add <- data.table::melt.data.table(add, 
                                       measure.vars = colnames(add)[-1], 
                                       variable.name = "col")
    # Inherit row clustering
    add[DT, c("rcl", "y"):= .(i.rcl, i.y), on= "row"]
    # Replace
    mat <- newdata
    DT <- add
  }

  #------------------------####
  # Clustering cols
  #------------------------####
  if(cluster_cols)
  {
    set.seed(3453)
    if(clustering_distance_cols %in% c("pearson", "spearman"))
    {
      .d <- as.dist(1 - cor(mat, use= "pairwise.complete.obs", method= clustering_distance_cols))
    }else
    {
      .d <- dist(t(mat), method = clustering_distance_cols)
    }
    ccl <- hclust(.d, method = clustering_method) # Hierarchical clustering
    cdend <- data.table::as.data.table(ggdendro::dendro_data(ccl, type = "rectangle", rotate= T)$segments) # Extract dendrogram
    cx <- data.table::data.table(ccl$labels, order(ccl$order))
    DT[cx, c("ccl", "x"):= .(1, i.V2), on= "col==V1"]
    if(!is.na(cutree_cols))
    {
      ccl <- data.table::data.table(names(cutree(ccl, cutree_cols)), cutree(ccl, cutree_cols))
      DT[ccl, ccl:= i.V2, on= "col==V1"]
    }
  }else
  {
    DT[, c("ccl", "x"):= .(1, .GRP), col] # No clustering
  }
  
  #------------------------####
  # Compute colors
  #------------------------####
  if(is.null(breaks))
    breaks <- seq(min(DT$value, na.rm= T), max(DT$value, na.rm= T), length.out = length(col))
  Cc <- circlize::colorRamp2(breaks, colors= col)
  DT[!is.na(value), Cc:= Cc(value)]
  DT[is.na(Cc), Cc:= "lightgrey"]

  #------------------------####
  # PLOT
  #------------------------####
  # Image
  im <- as.matrix(data.table::dcast(DT, y~x, value.var = "Cc"), 1)
  plot.new()
  rasterImage(im,
              xleft = 0,
              ybottom = 1,
              xright = 1,
              ytop = 0,
              interpolate = F)

  # Plot numbers
  if(display_numbers)
    text(DT$x/(max(DT$x))-(0.5/ncol(im)),
         DT$y/(max(DT$y))-(0.5/nrow(im)),
         display_numbers_FUN(if(is.null(newdata)) DT$value else DT$newdata),
         cex= display_numbers_cex,
         offset= 0)

  # Plot dendro and cuts
  if(cluster_rows & is.null(newdata))
    if(is.na(kmeans_k))
      segments(rdend$y/max(rdend$y)/10+1,
               1-(rdend$x/max(rdend$x)-(0.5/nrow(im))),
               rdend$yend/max(rdend$yend)/10+1,
               1-(rdend$xend/max(rdend$xend)-(0.5/nrow(im))),
               xpd= T)
    cr <- cumsum(unique(DT[, rcl, keyby= y])[, .N, rcl][, N])
    cr <- cr[-length(cr)]/nrow(im)
    if(length(cr)>0)
      segments(0, cr, 1, cr)

  if(cluster_cols)
    segments(cdend$x/max(cdend$x)-(0.5/ncol(im)),
             cdend$y/max(cdend$y)/10+1,
             cdend$xend/max(cdend$xend)-(0.5/ncol(im)),
             cdend$yend/max(cdend$yend)/10+1,
             xpd= T)
    cc <- cumsum(unique(DT[, ccl, keyby= x])[, .N, ccl][, N])
    cc <- cc[-length(cc)]/ncol(im)
    if(length(cc)>0)
      segments(cc, 0, cc, 1)

  # Plot axes
  if(show_colnames)
    axis(1,
         at= seq(ncol(im))/ncol(im)-(0.5/ncol(im)),
         labels = unique(DT[, col, keyby= x])$col,
         lty= 0,
         las= 2)
    
  if(show_rownames)
    axis(2,
         at= seq(nrow(im))/nrow(im)-(0.5/nrow(im)),
         labels = unique(DT[, row, keyby= y])$row,
         lty= 0,
         las= 2)

  # Plot legend
  xleft <- legend_pos[2]
  ybottom <- legend_pos[1]
  xright <- legend_pos[4]
  ytop <- legend_pos[3]
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
       ytop+par("usr")[4]*0.05,
       labels = legend_title,
       pos=4,
       xpd= T,
       cex= 0.8,
       offset= 0)
  
  invisible(DT)
}


