#' Heatmap
#'
#' This function plots a violin plot from either matrix,
#' data.table, data.frame or using formula.
#'
#' @param mat Matrix to plot. Row names and col names will be reported in return
#' @param newdata If specified, clusters using matrix and plots new data instead
#' @param breaks Heatmap breaks, which should be the same length as the color vector
#' @param show_rownames Show rownames? default= T
#' @param show_colnames Show colnames? default= T
#' @param col Color scale? default= c("cornflowerblue", "white", "red")
#' @param main Title, default is NA
#' @param auto_margins If T, optimal margins will be computed. Otherwise, takes existing ones
#' @param cluster_rows Should matrix rows be clustered? default is T
#' @param clustering_distance_rows Method for clustering distance rows. Can be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
#' @param cutree_rows Number of cuts for rows tree. default= 1z
#' @param show_cl_number_rows if cutree rows, should row cluster number be shown? default is TRUE
#' @param cluster_cols Should matrix cols be clustered? default is T
#' @param clustering_distance_cols Method for clustering distance cols. Can be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
#' @param cutree_cols Number of cuts for cols tree. default= 1z
#' @param show_cl_number_cols if cutree cols, should col cluster number be shown? default is TRUE
#' @param clustering_method Clustering method to use. Must be one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".
#' @param kmeans_k Number of row kmean clusters. If specified, takes over row clustering.
#' @param display_numbers Display numbers on heatmap? Default= F
#' @param display_numbers_FUN Function to apply before displaying numbers on heatmap.
#' @param display_numbers_cex cex display numbers
#' @param legend_title Character to plot as title legend
#' @param scale Scale matrix before clustering? Can be 'none' (default), 'column' or 'row'.
#' @param na_col Color for na values. Default= "lightgrey"
#' @param plot Should the heatmap be plotted? Default to T
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
                       main= NA,
                       auto_margins= T,
                       cluster_rows= T,
                       clustering_distance_rows= "euclidean",
                       cutree_rows = 1,
                       show_cl_number_rows= T,
                       cluster_cols = T,
                       clustering_distance_cols = "euclidean",
                       cutree_cols = 1,
                       show_cl_number_cols= T,
                       clustering_method = "complete",
                       kmeans_k = NA,
                       display_numbers= F,
                       display_numbers_FUN= function(x) round(x, 2),
                       display_numbers_cex= 0.5,
                       legend_title= "legend",
                       scale = "none",
                       na_col= "lightgrey",
                       plot= T)
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
  if(length(main)!=1)
    stop("length(main)!=1")

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
  if(cluster_rows & nrow(mat)>1)
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
  if(cluster_cols & ncol(mat)>1)
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
  DT[is.na(Cc), Cc:= na_col]
  
  #------------------------####
  # Compute numer plotting pos
  #------------------------####
  xpos <- 1/(max(DT$x, na.rm = T)*2) 
  xpos <- seq(xpos, 1-xpos, length.out= max(DT$x))
  ypos <- 1/(max(DT$y, na.rm = T)*2) 
  ypos <- seq(ypos, 1-ypos, length.out= max(DT$y))
  DT[, xplot:= xpos[x], x]
  DT[, yplot:= ypos[y], y]

  #------------------------####
  # PLOT
  #------------------------####
  if(plot)
  {
    # Margins
    if(auto_margins)
    {
      mBottom <- 0.5
      if(show_colnames)
        mBottom <- mBottom+max(strwidth(DT$col, "inches"))
      mLeft <- 0.5
      if(show_rownames)
        mLeft <- mLeft+max(strwidth(DT$row, "inches"))
      mTop <- 1
      if(cluster_cols)
        mTop <- mTop+0.5
      mRight <- strwidth(legend_title, "inches")+grconvertX(1, "lines", "inches")
      if(!is.na(kmeans_k))
        mRight <- mRight+grconvertX(1, "lines", "inches") else if(cluster_rows)
          mRight <- mRight+grconvertX(0.1, "npc", "inches")
      if(mRight<grconvertX(3, "lines", "inches"))
        mRight <- grconvertX(3, "lines", "inches")
      DT[1, mai:= .(c(mBottom, mLeft, mTop, mRight))]
      par(mai= c(mBottom, 
                 mLeft, 
                 mTop, 
                 mRight))
    }
    
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
      text(DT$xplot,
           DT$yplot,
           display_numbers_FUN(DT$value),
           cex= display_numbers_cex,
           offset= 0)
    
    # Plot dendro and cuts
    if(cluster_rows & is.null(newdata) & nrow(mat)>1)
      if(is.na(kmeans_k))
        segments(rdend$y/max(rdend$y)/10+1,
                 1-(rdend$x/max(rdend$x)-(0.5/nrow(im))),
                 rdend$yend/max(rdend$yend)/10+1,
                 1-(rdend$xend/max(rdend$xend)-(0.5/nrow(im))),
                 xpd= T)
    cr <- cumsum(unique(DT[, rcl, keyby= y])[, .N, rcl][, N])
    cr <- cr[-length(cr)]/nrow(im)
    if(length(cr)>0) # Plot cutree lines and cluster numbers
    {
      segments(0, cr, 1, cr)
      if(show_cl_number_rows)
        text(x= 1, 
             y= DT[, mean(y), rcl]$V1/(nrow(im))-(0.5/nrow(im)), 
             labels = DT[, mean(y), rcl]$rcl, 
             xpd= T, 
             pos= 4, 
             offset= 0.5,
             cex= 1.5)
    }
    
    if(cluster_cols & ncol(mat)>1)
      segments(cdend$x/max(cdend$x)-(0.5/ncol(im)),
               cdend$y/max(cdend$y)/10+1,
               cdend$xend/max(cdend$xend)-(0.5/ncol(im)),
               cdend$yend/max(cdend$yend)/10+1,
               xpd= T)
    cc <- cumsum(unique(DT[, ccl, keyby= x])[, .N, ccl][, N])
    cc <- cc[-length(cc)]/ncol(im)
    if(length(cc)>0) # Plot cutcol lines and cluster numbers
    {
      segments(cc, 0, cc, 1)
      if(show_cl_number_cols)
        text(x= DT[, mean(x), ccl]$V1/ncol(im)-0.5/ncol(im), 
             y= 1, 
             labels = DT[, mean(x), ccl]$ccl, 
             xpd= T, 
             pos= 3, 
             offset= 0.5,
             cex= 1.5)
    }
    
    # Plot axes
    if(show_colnames)
      axis(1,
           at= seq(ncol(im))/ncol(im)-(0.5/ncol(im)),
           labels = unique(DT[, col, keyby= x])$col,
           lty= 0,
           las= 2, 
           line= -1)
    
    if(show_rownames)
      axis(2,
           at= seq(nrow(im))/nrow(im)-(0.5/nrow(im)),
           labels = unique(DT[, row, keyby= y])$row,
           lty= 0,
           las= 2,
           line= -1)
    
    # Plot legend
    xleft <- grconvertX(1, "npc", "inches")+grconvertX(1, "lines", "inches")
    if(!is.na(kmeans_k))
      xleft <- xleft+grconvertX(1, "lines", "inches") else if(cluster_rows)
        xleft <- xleft+grconvertX(0.1, "npc", "inches")
    xright <- xleft+grconvertX(1, "lines", "inches")
    xleft <- grconvertX(xleft, "inches", "npc")
    xright <- grconvertX(xright, "inches", "npc")
    ybottom <- 0.7
    ytop <- grconvertY(1, "npc", "inches")-grconvertY(1, "chars", "inches")
    ytop <- grconvertY(ytop, "inches", "npc")
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
         1-0.5*strheight(legend_title),
         labels = legend_title,
         pos=4,
         xpd= T,
         cex= 0.8,
         offset= 0)
    
    # Title
    text(0.5, 
         y = ifelse(cluster_cols, 1.1, 1.05),
         main, 
         pos = 3,
         xpd= T)
  }
  
  invisible(DT)
}


