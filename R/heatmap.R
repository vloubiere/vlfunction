#' Heatmap
#'
#' This function plots a heatmap and returns an S3 object of type vl_heatmap_obj which can be modified and used for later plotting.
#' 
#' @param x Data to plot. Can be one of matrix, data.table, or formula.
#' @param cluster_rows Should rows be clustered? Default= T
#' @param kmeans_k Number of row kmean clusters. If specified, takes over row clustering.
#' @param clustering_distance_rows Method for clustering distance rows. Can be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
#' @param cutree_rows Number of cuts for rows tree. Default= 1L
#' @param cluster_cols Should columns be clustered? Default= T
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
#' @param display_numbers_FUN Function to apply before displaying numbers on heatmap.
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

#' @param x 
#' @param ... Extra parameters to be passed to ?vl_heatmap().
#' @describeIn vl_heatmap Just a wrapper that collapses DT to a suitable matrix
#' @export
vl_heatmap.data.table <- function(x, 
                                  rownames= NULL,
                                  rownames.value= NULL, 
                                  ...)
{
  mat <- data.table::as.matrix(x, 
                               rownames= rownames, 
                               rownames.value= rownames.value)
  vl_heatmap.matrix(mat, ...)
}

#' @describeIn vl_heatmap Default function using matrix as an input
#' @export
vl_heatmap.matrix <- function(x,
                              cluster_rows= T,
                              kmeans_k= NA,
                              clustering_distance_rows= "euclidean",
                              cutree_rows = 1,
                              cluster_cols = T,
                              clustering_distance_cols = "euclidean",
                              cutree_cols = 1,
                              clustering_method = "complete",
                              plot= T,
                              auto_margins= T,
                              breaks= seq(min(x, na.rm= T), max(x, na.rm= T), length.out= length(col)),
                              col= c("cornflowerblue", "white", "red"),
                              na_col= "lightgrey",
                              main= NA,
                              legend_title= NA,
                              show_rownames= T,
                              show_colnames= T,
                              show_row_clusters= T,
                              show_col_clusters= T,
                              show_row_dendrogram= T,
                              show_col_dendrogram= T,
                              show_legend= T,
                              display_numbers= F,
                              display_numbers_FUN= function(x) round(x, 2),
                              display_numbers_cex= 0.5)
{
  if(is.null(rownames(x)))
    rownames(x) <- seq(nrow(x))
  if(is.null(colnames(x)))
    colnames(x) <- seq(ncol(x))
  
  #------------------------####
  # Init informative result DT
  #------------------------####
  rows <- data.table(row= rownames(x), 
                     row_order= seq(nrow(x)),
                     row_cl= 1)
  cols <- data.table(col= colnames(x), 
                     col_order= seq(ncol(x)),
                     col_cl= 1)
  rcl <- NULL
  rdend <- NULL
  ccl <- NULL
  cdend <- NULL
  
  #------------------------####
  # Clustering rows
  #------------------------####
  if(cluster_rows)
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
                              rows[, row_order:= rcl$order]
                              # Cutree
                              rows[, row_cl:= cutree(rcl, cutree_rows)]
                              # Extract dend
                              dend <- ggdendro::dendro_data(rcl, 
                                                            type = "rectangle", 
                                                            rotate= T)
                              rdend <- data.table::as.data.table(dend$segments)
    }else
    {
      # Kmeans clustering
      rcl <- kmeans(x, centers = kmeans_k)
      rows[, row_order:= order(rcl$cluster)]
      rows[, row_cl:= rcl$cluster]
    }
  }
  
  #------------------------####
  # Clustering cols
  #------------------------####
  if(cluster_cols)
  {
    set.seed(3453)
    if(clustering_distance_cols %in% c("pearson", "spearman"))
      .d <- as.dist(1 - cor(x, 
                            use= "pairwise.complete.obs", 
                            method= clustering_distance_cols)) else
                              .d <- dist(t(x), method = clustering_distance_cols)
                            # Hierarchical clustering
                            ccl <- hclust(.d, method = clustering_method)
                            cols[, col_order:= ccl$order]
                            # Cutree
                            cols[, col_cl:= cutree(ccl, cutree_cols)]
                            # Extract dend
                            dend <- ggdendro::dendro_data(ccl,
                                                          type = "rectangle", 
                                                          rotate= T)
                            cdend <- data.table::as.data.table(dend$segments)
  }
  
  #------------------------####
  # PLOT
  #------------------------####
  obj <- list(x= x,
              rows= rows,
              cols= cols,
              rcl= rcl,
              ccl= ccl,
              rdend= rdend,
              cdend= cdend)
  setattr(obj, "class", c("vl_heatmap", "list"))
  if(plot)
  {
    pl <- match.call()
    pl$obj <- obj
    pl$x <- NULL
    pl$kmeans_k <- pl$clustering_distance_rows <- NULL
    pl$clustering_distance_cols <- NULL
    pl$clustering_method <- NULL
    pl[[1]] <- quote(plot.vl_heatmap)
    eval(pl)
  }
  
  invisible(obj)
}


#' Title
#'
#' @param obj vl_heatmap object containing the parent clustering
#' @param add matrix that will inhiherit obj clustering
#' @param inherit_row_cl Should row clustering be inherited? Default is TRUE if rownames match
#' @param inherit_col_cl Should col clustering be inherited? Default is TRUE if colnames match
#' @param ... Extra arguments to be passed to vl_heatmap
#' @return a vl_heatmap object whose clustering is inherited from obj
#' @describeIn vl_heatmap Inherit clustering to another matrix
#' @export
vl_heatmap.vl_heatmap <- function(obj, 
                                  add,
                                  inherit_row_cl= T,
                                  inherit_col_cl= T, 
                                  ...)
{
  cl <- vl_heatmap(add, plot= F, ...)
  
  if(inherit_row_cl & identical(rownames(add), rownames(obj$x)))
  {
    cl$rows <- obj$rows
    cl$rcl <- obj$rcl
    cl$rdend <- obj$rdend
  }else
    inherit_row_cl <- F
  if(inherit_col_cl & identical(colnames(add), colnames(obj$x)))
  {
    cl$cols <- obj$cols
    cl$ccl <- obj$ccl
    cl$cdend <- obj$cdend
  }else
    inherit_col_cl <- F
  if(!inherit_row_cl & !inherit_col_cl)
    stop("add should share either rownames or colnames with the original matrix from obj")
  
  #SAVE
  invisible(cl)
}

# Default plotting function
#' @describeIn vl_heatmap Default plotting function
#' @export
plot.vl_heatmap <- function(obj,
                            cluster_rows= T,
                            cluster_cols= T,
                            cutree_rows = 1,
                            cutree_cols = 1,
                            auto_margins= T,
                            breaks= NULL,
                            col= c("cornflowerblue", "white", "red"),
                            na_col= "lightgrey",
                            main= NA,
                            show_legend= T,
                            legend_title= NA,
                            show_rownames= T,
                            show_colnames= T,
                            show_row_clusters= T,
                            show_col_clusters= T,
                            show_row_dendrogram= T,
                            show_col_dendrogram= T,
                            display_numbers= F,
                            display_numbers_FUN= function(x) round(x, 2),
                            display_numbers_cex= 0.5)
{
  list2env(obj, environment())
  # Checks
  if(!cluster_rows | is.null(rcl))
  {
    cluster_rows <- F
    show_row_clusters <- F
    show_row_dendrogram <- F
  }else if(class(rcl)=="kmeans")
    show_row_dendrogram <- F
  if(!cluster_cols | is.null(ccl))
  {
    cluster_cols <- F
    show_col_clusters <- F
    show_col_dendrogram <- F
  }
  
  # Margins
  if(auto_margins)
  {
    bot <- 1
    left <- 1
    top <- 1.5
    right <- 1
    if(show_colnames)
      bot <- bot+grconvertY(max(strwidth(cols$col, "inches")), "inches", "lines")
    if(show_rownames)
      left <- left+grconvertX(max(strwidth(rows$row, "inches")), "inches", "lines")
    if(is.character(main))
      top <- top+grconvertY(strheight(main, cex = 0.8, units = "inches"), "inches", "lines")
    if(show_col_dendrogram)
      top <- top+1.5 else if(show_col_clusters)
          top <- top+1
    if(show_row_dendrogram)
      right <- right+1.5 else if(show_row_clusters)
        right <- right+1
    leg_width <- strwidth(legend_title, units = "inches", cex= 0.8)+par("cin")[1]
    leg_width <- grconvertX(leg_width, "inches", "lines")
    if(leg_width>3.5)
      right <- right+leg_width else if(show_legend)
        right <- right+3.5
    par(mar= c(bot, left, top, right))
  }
  
  #----------------------------------#
  # Heatmap
  #----------------------------------#
  # Order matrix
  if(cluster_rows)
    x <- x[(rows$row_order),]
  if(cluster_cols)
    x <- x[,(cols$col_order)]
  # Breaks
  if(is.null(breaks))
    breaks <- seq(min(x, na.rm= T), max(x, na.rm= T), length.out= length(col))
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
    text(c(col(x)),
         rev(c(row(x))),
         display_numbers_FUN(c(x)),
         cex= display_numbers_cex,
         offset= 0)
  
  # Cluser lines
  if(show_row_clusters)
  {
    if(class(rcl)=="hclust" && cutree_rows>max(rows$row_cl))
      rows[, row_cl:= cutree(rcl, cutree_rows)]
    pos <- cumsum(rev(rows[(row_order), .N, row_cl]$N))
    abline(h= pos[-length(pos)]+0.5)
    text(x = par("usr")[2], 
         y= pos-diff(c(0, pos))/2, 
         labels = rev(unique(rows[(row_order), row_cl])),
         pos= 4,
         xpd= T)
  }
  if(show_col_clusters)
  {
    if(cutree_cols>max(cols$col_cl))
      cols[, col_cl:= cutree(ccl, cutree_cols)]
    pos <- cumsum(cols[(col_order), .N, col_cl]$N)
    abline(v= pos[-length(pos)]+0.5)
    text(x = pos-diff(c(0, pos))/2,
         y= par("usr")[4],
         labels = unique(cols[(col_order), col_cl]),
         pos= 3,
         xpd= T)
  }

  #----------------------------------#
  # Margins legends
  #----------------------------------#
  # Margin lines width and height in user coordinates
  mar.lw <- diff(grconvertX(c(0,1), "lines", "user"))
  mar.lh <- diff(grconvertY(c(0,1), "lines", "user"))

  # Title
  if(!is.na(main))
    mtext(main,
          line = ifelse(show_col_dendrogram, 1.75, 0.25))

  # Plot dendro and cuts
  if(show_row_dendrogram)
    segments(rdend$y/diff(range(rdend[,c(y, yend)]))*mar.lw*1.5+par("usr")[2],
             par("usr")[4]-rdend$x+0.5,
             rdend$yend/diff(range(rdend[,c(y, yend)]))*mar.lw*1.5+par("usr")[2],
             par("usr")[4]-rdend$xend+0.5,
             xpd= T)
  if(show_col_dendrogram)
    segments(cdend$x,
             cdend$y/diff(range(cdend[,c(y, yend)]))*mar.lh*1.5+par("usr")[4],
             cdend$xend,
             cdend$yend/diff(range(cdend[,c(y, yend)]))*mar.lh*1.5+par("usr")[4],
             xpd= T)

  # Plot axes
  if(show_colnames)
    axis(1,
         at= seq(ncol(x)),
         labels = colnames(x),
         lwd= 0,
         las= 2,
         line= -0.5)
  if(show_rownames)
    axis(2,
         at= seq(nrow(x)),
         labels = rev(rownames(x)),
         lwd= 0,
         las= 2,
         line= -0.5)

  # Plot legend
  if(show_legend)
  {
    left <- par("usr")[2]+mar.lw
    if(show_row_dendrogram)
      left <- left+mar.lw*1.5
    top <- par("usr")[4]-mar.lh*2
    width <- mar.lw
    vl_heatkey(breaks,
               col,
               left,
               top,
               mar.lh*6,
               mar.lw,
               main= legend_title)
  }
}
