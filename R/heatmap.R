#' Heatmap
#'
#' This function plots a heatmap and returns an S3 object of type vl_heatmap_obj which can be modified and used for later plotting.
#' 
#' @param x Data to plot. Can be one of matrix, data.table, or formula.
#' @param cluster_rows Should rows be clustered? Default= T
#' @param row_clusters Vector of row clusters (overwritten by clustering, Default= 1)
#' @param row_clusters_col Row clusters colors (default is to use grey palette)
#' @param kmeans_k Number of row kmean clusters. If specified, takes over row clustering.
#' @param clustering_distance_rows Method for clustering distance rows. Can be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
#' @param cutree_rows Number of cuts for rows tree. Default= 1L
#' @param cluster_cols Should columns be clustered? Default= T
#' @param col_clusters Vector of col clusters (overwritten by clustering, Default= 1)
#' @param col_clusters_col Col clusters colors (default is to use grey palette)
#' @param clustering_distance_cols Method for clustering distance cols. Can be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
#' @param cutree_cols  Number of cuts for cols tree. default= 1L
#' @param clustering_method Clustering method to use. Must be one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".
#' @param plot Should the heatmap be plotted? Default= T
#' @param auto_margins Should auto margins be computed and used?
#' @param breaks Heatmap breaks, which should be the same length as the color vector
#' @param col Color scale. default= c("cornflowerblue", "white", "red")
#' @param na_col Color for na values. Default= "lightgrey"
#' @param main Title. Default= NA
#' @param legend_title Character to plot as title legend
#' @param show_rownames Plot row names? Default= T
#' @param show_colnames Plot col names? Default= T
#' @param show_row_clusters Plot row cluster names/lines?
#' @param show_col_clusters Plot col cluster names/lines?
#' @param show_row_dendrogram Plot row dendrogram?
#' @param show_col_dendrogram Plot col dendrogram?
#' @param show_legend Should the legend be plotted?
#' @param display_numbers Display numbers on heatmap? Default= F
#' @param display_numbers_matrix Matrix of numbers to be displayed. useful for fine tuning. Default= x
#' @param display_numbers_cex cex display numbers 
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
#' vl_heatmap(test, cutree_rows = 3, cutree_cols = 2, show_col_clusters = F, legend_title = "log2FoldChange", main= "test")
#' 
#' @return An object of class vl_heatmap_obj containing clustering data and can be used for plotting using plot(obj, ...)
#' @export
vl_heatmap <- function(x, ...) UseMethod("vl_heatmap")

#' @describeIn vl_heatmap Just a wrapper that collapses DT to a suitable matrix
#' @export
vl_heatmap.data.table <- function(x, 
                                  rownames= NULL,
                                  rownames.value= NULL, 
                                  ...)
{
  mat <- as.matrix(x, 
                   rownames= rownames, 
                   rownames.value= rownames.value)
  vl_heatmap.matrix(mat, ...)
}

#' @describeIn vl_heatmap Default function using matrix as an input
#' @export
vl_heatmap.matrix <- function(x,
                              cluster_rows= TRUE,
                              row_clusters= 1,
                              row_clusters_col= NULL,
                              kmeans_k= NA,
                              clustering_distance_rows= "euclidean",
                              cutree_rows = 1,
                              cluster_cols = TRUE,
                              col_clusters= 1,
                              col_clusters_col= NULL,
                              clustering_distance_cols = "euclidean",
                              cutree_cols = 1,
                              clustering_method = "complete",
                              plot= TRUE,
                              auto_margins= FALSE,
                              breaks= seq(min(x, na.rm= T), max(x, na.rm= TRUE), length.out= length(col)),
                              col= c("cornflowerblue", "white", "red"),
                              na_col= "lightgrey",
                              main= NA,
                              legend_title= NA,
                              show_rownames= TRUE,
                              show_colnames= TRUE,
                              show_row_clusters= length(unique(rows$cl))>1,
                              show_col_clusters= length(unique(cols$cl))>1,
                              show_row_dendrogram= TRUE,
                              show_col_dendrogram= TRUE,
                              show_legend= TRUE,
                              display_numbers= FALSE,
                              display_numbers_matrix,
                              display_numbers_cex= 1)
{
  if(is.null(rownames(x)))
    rownames(x) <- seq(nrow(x))
  if(is.null(colnames(x)))
    colnames(x) <- seq(ncol(x))
  if(!is.factor(row_clusters))
    row_clusters <- factor(row_clusters)
  if(!is.factor(col_clusters))
    col_clusters <- factor(col_clusters)
  if(is.null(row_clusters_col))
    if(cluster_rows && length(levels(row_clusters))==1)
      row_clusters_col <- grDevices::gray.colors(cutree_rows) else
        row_clusters_col <- grDevices::gray.colors(length(unique(row_clusters)))
  if(is.null(col_clusters_col))
    if(cluster_cols && length(levels(col_clusters))==1)
      col_clusters_col <- grDevices::gray.colors(cutree_cols) else
        col_clusters_col <- grDevices::gray.colors(length(unique(col_clusters)))
        
  
  #------------------------####
  # Init informative result DT
  #------------------------####
  rows <- data.table(name= rownames(x),
                     cl= row_clusters)
  rows[, order:= order(cl)]
  cols <- data.table(name= colnames(x),
                     cl= col_clusters)
  cols[, order:= order(cl)]
  
  rcl <- NULL
  rdend <- NULL
  ccl <- NULL
  cdend <- NULL
  
  #------------------------####
  # Clustering rows
  #------------------------####
  if(cluster_rows && nrow(x)>1 && length(levels(row_clusters))==1)
  {
    set.seed(3453)
    if(is.na(kmeans_k))
    {
      if(clustering_distance_rows %in% c("pearson", "spearman"))
        .d <- as.dist(1 - cor(t(x), 
                              use= "pairwise.complete.obs", 
                              method= clustering_distance_rows)) else
                                .d <- dist(x, method = clustering_distance_rows)
                              # Hierarchical clustering
                              rcl <- hclust(.d, method = clustering_method) 
                              rows[, order:= rcl$order]
                              # Cutree
                              rows[, cl:= cutree(rcl, cutree_rows)]
                              # Extract dend
                              dend <- ggdendro::dendro_data(rcl, 
                                                            type = "rectangle", 
                                                            rotate= T)
                              rdend <- data.table::as.data.table(dend$segments)
    }else
    {
      # Kmeans clustering
      rcl <- kmeans(x, centers = kmeans_k)
      rows[, order:= order(rcl$cluster)]
      rows[, cl:= rcl$cluster]
    }
  }
  # Add cluster color and y pos
  rows[, col:= row_clusters_col[cl]]
  rows[(order), y:= rev(.I)]
  
  #------------------------####
  # Clustering cols
  #------------------------####
  if(cluster_cols && ncol(x)>1 && length(levels(col_clusters))==1)
  {
    set.seed(3453)
    if(clustering_distance_cols %in% c("pearson", "spearman"))
      .d <- as.dist(1 - cor(x, 
                            use= "pairwise.complete.obs", 
                            method= clustering_distance_cols)) else
                              .d <- dist(t(x), method = clustering_distance_cols)
                            # Hierarchical clustering
                            ccl <- hclust(.d, method = clustering_method)
                            cols[, order:= ccl$order]
                            # Cutree
                            cols[, cl:= cutree(ccl, cutree_cols)]
                            # Extract dend
                            dend <- ggdendro::dendro_data(ccl,
                                                          type = "rectangle", 
                                                          rotate= T)
                            cdend <- data.table::as.data.table(dend$segments)
  }
  # Add cluster color and x pos
  cols[, col:= col_clusters_col[cl]]
  cols[(order), x:= .I]
  
  #------------------------####
  # PLOT
  #------------------------####
  obj <- mget(ls())
  setattr(obj, "class", c("vl_heatmap", "list"))
  if(plot)
    plot.vl_heatmap(obj)
  
  # RETURN
  invisible(obj)
}

#' @export
plot.vl_heatmap <- function(obj)
{
  list2env(obj, environment())
  
  # Margins
  if(auto_margins)
  {
    bot <- 1
    left <- 1
    if(show_colnames)
      bot <- bot+grconvertY(max(strwidth(cols$name, "inches")), "inches", "lines")
    if(show_rownames)
      left <- left+grconvertX(max(strwidth(rows$name, "inches")), "inches", "lines")
    par(mar= c(bot, left, 5, 7))
  }
  
  #----------------------------------#
  # Heatmap
  #----------------------------------#
  # Order matrix
  x <- x[(rows$order), , drop=F]
  x <- x[, (cols$order), drop=F]
    
  # Palette
  Cc <- circlize::colorRamp2(breaks, 
                             colors= col)
  # Image
  im <- x
  im[!is.na(im)] <- Cc(im[!is.na(im)])
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
  box(lwd= 0.25)
  
  # Plot numbers
  if(display_numbers)
  {
    if(missing(display_numbers_matrix))
      display_numbers_matrix <- x else if(identical(dim(x), dim(display_numbers_matrix)))
      {
        display_numbers_matrix <- display_numbers_matrix[(rows$order), , drop=F]
        display_numbers_matrix <- display_numbers_matrix[, (cols$order), drop=F]
      }else
        stop("dim(x) and dim(display_numbers_matrix) should be identical!")
    
    text(c(col(display_numbers_matrix)),
         rev(c(row(display_numbers_matrix))),
         c(display_numbers_matrix),
         cex= display_numbers_cex,
         offset= 0) 
  }

  #----------------------------------#
  # Margins legends
  #----------------------------------#
  # Margin lines width and height in user coordinates
  mar.lw <- diff(grconvertX(c(0,1), "lines", "user"))
  mar.lh <- diff(grconvertY(c(0,1), "lines", "user"))

  # Cluser lines
  if(show_row_clusters)
  {
    pos <- rows[rev(order)][, .(y1= .N), .(cl, col)]
    pos[, y1:= cumsum(y1)+0.5]
    pos[, y0:= data.table::shift(y1, 1, fill = 0.5)]
    abline(h= data.table::first(pos$y1, nrow(pos)-1))
    rleft <- par("usr")[2]+mar.lw/10
    rect(rleft,
         pos$y0,
         rleft+mar.lw,
         pos$y1,
         col= pos$col,
         xpd= T,
         border= NA)
    text(rleft+mar.lw/2,
         rowMeans(pos[, .(y0, y1)]),
         pos$cl,
         srt= 270,
         xpd= T,
         adj= 0.5)
  }
  if(show_col_clusters)
  {
    pos <- cols[(order)][, .(x1= .N), .(cl, col)]
    pos[, x1:= cumsum(x1)+0.5]
    pos[, x0:= data.table::shift(x1, 1, fill = 0.5)]
    abline(v= data.table::first(pos$x1, nrow(pos)-1))
    ctop <- par("usr")[4]+mar.lh/10
    rect(pos$x0,
         ctop,
         pos$x1,
         ctop+mar.lh,
         col= pos$col,
         xpd= T,
         border= NA)
    text(rowMeans(pos[, .(x0, x1)]),
         ctop+mar.lh/2,
         pos$cl,
         xpd= T,
         adj= 0.5)
  }

  # Title
  if(!is.na(main))
    title(main= main)

  # Plot dendrograms
  if(show_row_dendrogram & !is.null(rdend))
  {
    if(show_row_clusters)
      d.left <- rleft+mar.lw else
        d.left <- par("usr")[2]
    segments(rdend$y/diff(range(rdend[,c(y, yend)]))*mar.lw*1.5+d.left,
             par("usr")[4]-rdend$x+0.5,
             rdend$yend/diff(range(rdend[,c(y, yend)]))*mar.lw*1.5+d.left,
             par("usr")[4]-rdend$xend+0.5,
             xpd= T)
  }
  if(show_col_dendrogram & !is.null(cdend))
  {
    if(show_col_clusters)
      d.bot <- ctop+mar.lh else
        d.bot <- par("usr")[4]
    segments(cdend$x,
             cdend$y/diff(range(cdend[,c(y, yend)]))*mar.lh*1.5+d.bot,
             cdend$xend,
             cdend$yend/diff(range(cdend[,c(y, yend)]))*mar.lh*1.5+d.bot,
             xpd= T)
  }

  # Plot axes
  if(show_colnames)
    axis(1,
         at= seq(ncol(x)),
         labels = colnames(x),
         lwd= 0)
  if(show_rownames)
    axis(2,
         at= seq(nrow(x)),
         labels = rev(rownames(x)),
         lwd= 0)

  # Plot legend
  if(show_legend)
  {
    left <- par("usr")[2]+mar.lw
    if(show_row_clusters)
      left <- left+mar.lw*1
    if(show_row_dendrogram && !is.null(rdend))
      left <- left+mar.lw*1.5
    top <- par("usr")[4]-mar.lh*2
    width <- mar.lw
    vl_heatkey(breaks,
               main.cex= 0.8,
               col,
               left,
               top,
               mar.lh*6,
               mar.lw,
               main= legend_title)
  }
}
