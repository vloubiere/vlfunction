#' Heatmap
#'
#' This function plots a heatmap and returns an S3 object of type vl_heatmap_obj which can be modified and used for later plotting.
#' 
#' @param x Data to plot. Can be one of matrix, data.table, or formula.
#' @param cluster.rows Should rows be clustered? Default= T
#' @param row.clusters Vector of row clusters (overwritten by clustering, Default= 1)
#' @param row.clusters.col Row clusters colors (default is to use grey palette)
#' @param kmeans.k Number of row kmean clusters. If specified, takes over row clustering.
#' @param clustering.distance.rows Method for clustering distance rows. Can be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
#' @param cutree.rows Number of cuts for rows tree. Default= 1L
#' @param cluster.cols Should columns be clustered? Default= T
#' @param col.clusters Vector of col clusters (overwritten by clustering, Default= 1)
#' @param col.clusters.col Col clusters colors (default is to use grey palette)
#' @param clustering.distance.cols Method for clustering distance cols. Can be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
#' @param cutree.cols  Number of cuts for cols tree. default= 1L
#' @param clustering.method Clustering method to use. Must be one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".
#' @param plot Should the heatmap be plotted? Default= T
#' @param breaks Heatmap breaks, which should be the same length as the color vector
#' @param col Color scale. default= c("cornflowerblue", "white", "red")
#' @param na.col Color for na values. Default= "lightgrey"
#' @param grid Should grid be drawn? Default= F
#' @param main Title. Default= NA
#' @param legend.title Character to plot as title legend
#' @param legend.cex cex expansion factor applying to the whole legend
#' @param show.rownames Plot row names? Default= T
#' @param show.colnames Plot col names? Default= T
#' @param tilt.colnames Should colnames be tilted? default= F
#' @param show.row.clusters Plot row cluster names/lines?
#' @param row.clusters.pos One of "left" or "right"
#' @param row.cluster.line.col Color of the line separting row clusters
#' @param show.col.clusters Plot col cluster names/lines?
#' @param show.row.dendrogram Plot row dendrogram?
#' @param show.col.dendrogram Plot col dendrogram?
#' @param show.legend Should the legend be plotted?
#' @param display.numbers Display numbers on heatmap? Default= F
#' @param display.numbers.matrix Matrix of numbers to be displayed. useful for fine tuning. Default= x
#' @param display.numbers.cex cex display numbers 
#' @param box.lwd Line width of the box around the heatmap
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
#' res <- vl_heatmap(test,
#'                   kmeans.k = 4,
#'                   cutree.rows = 3,
#'                   cutree.cols = 5,
#'                   legend.title = "log2FoldChange",
#'                   main= "test")
#' vl_heatmap(test,
#'            cutree.rows = 3,
#'            cutree.cols = 2,
#'            show.col.clusters = F,
#'            legend.title = "log2FoldChange",
#'            main= "test")
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
                              cluster.rows= TRUE,
                              row.clusters= 1,
                              row.clusters.col= NULL,
                              kmeans.k= NA,
                              clustering.distance.rows= "euclidean",
                              cutree.rows = 1,
                              cluster.cols = TRUE,
                              col.clusters= 1,
                              col.clusters.col= NULL,
                              clustering.distance.cols = "euclidean",
                              cutree.cols = 1,
                              clustering.method = "complete",
                              plot= TRUE,
                              breaks= seq(min(x, na.rm= T), max(x, na.rm= TRUE), length.out= length(col)),
                              col= c("cornflowerblue", "white", "red"),
                              na.col= "lightgrey",
                              grid= FALSE,
                              main= NA,
                              legend.title= NA,
                              legend.cex= 1,
                              show.rownames= TRUE,
                              show.colnames= TRUE,
                              tilt.colnames= FALSE,
                              show.row.clusters= length(unique(rows$cl))>1,
                              row.clusters.pos= "right",
                              row.cluster.line.col= "black",
                              show.col.clusters= length(unique(cols$cl))>1,
                              show.row.dendrogram= TRUE,
                              show.col.dendrogram= TRUE,
                              show.legend= TRUE,
                              display.numbers= FALSE,
                              display.numbers.matrix,
                              display.numbers.cex= 1,
                              box.lwd= .25)
{
  if(is.null(rownames(x)))
    rownames(x) <- seq(nrow(x))
  if(is.null(colnames(x)))
    colnames(x) <- seq(ncol(x))
  if(!is.factor(row.clusters))
    row.clusters <- factor(row.clusters)
  if(!is.factor(col.clusters))
    col.clusters <- factor(col.clusters)
  if(is.null(row.clusters.col))
  {
    if(cluster.rows && !is.na(kmeans.k))
    {
      row.clusters.col <- grDevices::gray.colors(kmeans.k)
    }else if(cluster.rows && length(levels(row.clusters))==1)
    {
      row.clusters.col <- grDevices::gray.colors(cutree.rows)
    }else
      row.clusters.col <- grDevices::gray.colors(length(unique(row.clusters)))
  }
  if(is.null(col.clusters.col))
  {
    if(cluster.cols && length(levels(col.clusters))==1)
      col.clusters.col <- grDevices::gray.colors(cutree.cols) else
        col.clusters.col <- grDevices::gray.colors(length(unique(col.clusters)))
  }
  if(row.clusters.pos=="right")
    show.rownames <- F
    
  # Init informative result DT ----
  rows <- data.table(name= rownames(x),
                     cl= row.clusters)
  rows[, order:= order(cl)]
  cols <- data.table(name= colnames(x),
                     cl= col.clusters)
  cols[, order:= order(cl)]
  
  rcl <- NULL
  rdend <- NULL
  ccl <- NULL
  cdend <- NULL
  
  # Clustering rows ----
  if(cluster.rows && nrow(x)>1 && length(levels(row.clusters))==1)
  {
    set.seed(3453)
    if(is.na(kmeans.k))
    {
      if(clustering.distance.rows %in% c("pearson", "spearman"))
        .d <- as.dist(1 - cor(t(x), 
                              use= "pairwise.complete.obs", 
                              method= clustering.distance.rows)) else
                                .d <- dist(x, method = clustering.distance.rows)
                              # Hierarchical clustering
                              rcl <- hclust(.d, method = clustering.method) 
                              rows[, order:= rcl$order]
                              # Cutree
                              rows[, cl:= cutree(rcl, cutree.rows)]
                              # Extract dend
                              dend <- ggdendro::dendro_data(rcl, 
                                                            type = "rectangle", 
                                                            rotate= T)
                              rdend <- data.table::as.data.table(dend$segments)
    }else
    {
      # Kmeans clustering
      rcl <- kmeans(x, centers = kmeans.k)
      rows[, order:= order(rcl$cluster)]
      rows[, cl:= rcl$cluster]
    }
  }
  # Add cluster color and y pos
  rows[, col:= row.clusters.col[cl]]
  rows[(order), y:= rev(.I)]
  
  # Clustering cols ----
  if(cluster.cols && ncol(x)>1 && length(levels(col.clusters))==1)
  {
    set.seed(3453)
    if(clustering.distance.cols %in% c("pearson", "spearman"))
      .d <- as.dist(1 - cor(x, 
                            use= "pairwise.complete.obs", 
                            method= clustering.distance.cols)) else
                              .d <- dist(t(x), method = clustering.distance.cols)
                            # Hierarchical clustering
                            ccl <- hclust(.d, method = clustering.method)
                            cols[, order:= ccl$order]
                            # Cutree
                            cols[, cl:= cutree(ccl, cutree.cols)]
                            # Extract dend
                            dend <- ggdendro::dendro_data(ccl,
                                                          type = "rectangle", 
                                                          rotate= T)
                            cdend <- data.table::as.data.table(dend$segments)
  }
  # Add cluster color and x pos
  cols[, col:= col.clusters.col[cl]]
  cols[(order), x:= .I]
  
  # PLOT ----
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
  
  # Heatmap ----
  # Order matrix
  x <- x[(rows$order), , drop=F]
  x <- x[, (cols$order), drop=F]
  # Palette
  Cc <- circlize::colorRamp2(breaks, 
                             colors= col)
  # Image
  im <- x
  im[!is.na(im)] <- Cc(im[!is.na(im)])
  im[is.na(im)] <- na.col
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
  if(box.lwd)
    box(lwd= box.lwd)
  
  # plot grid ----
  if(grid)
  {
    segments(c(0.5, seq(ncol(x))+0.5),
             0.5,
             c(0.5, seq(ncol(x))+0.5),
             nrow(x)+0.5,
             xpd= T)
    segments(0.5,
             c(0.5, seq(nrow(x))+0.5),
             ncol(x)+0.5,
             c(0.5, seq(nrow(x))+0.5),
             xpd= T)
  }
  # Plot numbers ----
  if(display.numbers)
  {
    if(missing(display.numbers.matrix))
      display.numbers.matrix <- x else if(identical(dim(x), dim(display.numbers.matrix)))
      {
        display.numbers.matrix <- display.numbers.matrix[(rows$order), , drop=F]
        display.numbers.matrix <- display.numbers.matrix[, (cols$order), drop=F]
      }else
        stop("dim(x) and dim(display.numbers.matrix) should be identical!")
    
    text(c(col(display.numbers.matrix)),
         rev(c(row(display.numbers.matrix))),
         c(display.numbers.matrix),
         cex= display.numbers.cex,
         offset= 0) 
  }

  # Margins legends ----
  # Margin lines width and height in user coordinates
  mar.lw <- diff(grconvertX(c(0,1), "lines", "user"))
  mar.lh <- diff(grconvertY(c(0,1), "lines", "user"))

  # Cluster lines ----
  if(show.row.clusters)
  {
    pos <- rows[rev(order)][, .(y1= .N), .(cl, col)]
    pos[, y1:= cumsum(y1)+0.5]
    pos[, y0:= data.table::shift(y1, 1, fill = 0.5)]
    abline(h= data.table::first(pos$y1, nrow(pos)-1),
           col= row.cluster.line.col)
    if(row.clusters.pos=="right")
    {
      rleft <- par("usr")[2]+mar.lw/10
      rect(rleft,
           pos$y0,
           rleft+mar.lw,
           pos$y1,
           col= pos$col,
           xpd= NA,
           border= NA)
      text(rleft+mar.lw/2,
           rowMeans(pos[, .(y0, y1)]),
           pos$cl,
           srt= 270,
           xpd= NA,
           adj= 0.5)
    }else if(row.clusters.pos=="left")
    {
      text(par("usr")[1]-diff(grconvertX(c(0, par("mgp")[1]+.5), "line", "user")),
           rowMeans(pos[, .(y0, y1)]),
           offset= 0,
           pos$cl,
           xpd= NA,
           adj= 1,
           cex= par("cex.lab"))
    }
  }
  if(show.col.clusters)
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
         xpd= NA,
         border= NA)
    text(rowMeans(pos[, .(x0, x1)]),
         ctop+mar.lh/2,
         pos$cl,
         xpd= NA,
         adj= 0.5)
  }

  # Title ----
  if(!is.na(main))
    title(main= main)

  # Plot dendrograms ----
  if(show.row.dendrogram & !is.null(rdend))
  {
    if(show.row.clusters)
      d.left <- rleft+mar.lw else
        d.left <- par("usr")[2]
    segments(rdend$y/diff(range(rdend[,c(y, yend)]))*mar.lw*1.5+d.left,
             par("usr")[4]-rdend$x+0.5,
             rdend$yend/diff(range(rdend[,c(y, yend)]))*mar.lw*1.5+d.left,
             par("usr")[4]-rdend$xend+0.5,
             xpd= NA)
  }
  if(show.col.dendrogram & !is.null(cdend))
  {
    if(show.col.clusters)
      d.bot <- ctop+mar.lh else
        d.bot <- par("usr")[4]
    segments(cdend$x,
             cdend$y/diff(range(cdend[,c(y, yend)]))*mar.lh*1.5+d.bot,
             cdend$xend,
             cdend$yend/diff(range(cdend[,c(y, yend)]))*mar.lh*1.5+d.bot,
             xpd= NA)
  }

  # Plot axes ----
  if(show.colnames)
  {
    if(tilt.colnames)
      vl_tilt_xaxis(x = seq(ncol(x)),
                    y = grconvertY(grconvertY(0, "npc", "inch")-grconvertY(par("mgp")[2], "line", "inch"), "inch", "user"),
                    labels = colnames(x),
                    offset= 0.25*par("mgp")[2],
                    pos= 2,
                    xpd= NA,
                    cex= par("cex.axis")) else
                      axis(1,
                           at= seq(ncol(x)),
                           labels = colnames(x),
                           lwd= 0)
  }
  if(show.rownames)
    axis(2,
         at= seq(nrow(x)),
         labels = rev(rownames(x)),
         lwd= 0)

  # Plot legend ----
  if(show.legend)
  {
    left <- par("usr")[2]+mar.lw
    if(show.row.clusters)
      left <- left+mar.lw*1
    if(show.row.dendrogram && !is.null(rdend))
      left <- left+mar.lw*1.5
    top <- par("usr")[4]-mar.lh
    vl_heatkey(breaks = breaks,
               col= col,
               left= left,
               top= top,
               height= mar.lh*6*legend.cex,
               width= mar.lw*legend.cex,
               ticks.cex= 0.6*legend.cex,
               main.cex= 0.8*legend.cex,
               main= legend.title)
  }
}