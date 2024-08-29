#' Heatmap
#'
#' This function plots a heatmap and returns an S3 object of type vl_heatmap_obj which can be modified and used for later plotting.
#' 
#' @param x Data to plot. Can be one of matrix, data.table, or formula.
#' @param cluster.rows Should rows be clustered? Default= TRUE.
#' @param clustering.distance.rows Method for clustering distance rows. Can be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall".
#' @param cutree.rows Number of cuts for rows tree. Default= 1L.
#' @param row.clusters Vector of row clusters (overwritten by clustering). Default= 1L (no clusters).
#' @param kmeans.k Number of row kmean clusters. If specified, takes over row clustering. Default= NA.
#' @param cluster.cols Should columns be clustered? Default= TRUE.
#' @param clustering.distance.cols Method for clustering distance cols. Can be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall".
#' @param cutree.cols  Number of cuts for cols tree. Default= 1L.
#' @param col.clusters Vector of col clusters (overwritten by clustering). Default= 1L (no clusters).
#' @param clustering.method Clustering method to use. Must be one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".
#' @param plot Should the heatmap be plotted? Default= T.
#' @param breaks Heatmap breaks, which should be the same length as the color vector.
#' @param col Color scale. default= c("cornflowerblue", "white", "red").
#' @param na.col Color for na values. Default= "lightgrey".
#' @param useRaster Should rasterImage be used? Default= TRUE.
#' @param grid Should grid be drawn? Default= FALSE.
#' @param main Title. Default= NA.
#' @param show.rownames Plot row names? Default= TRUE.
#' @param show.colnames Plot col names? Default= TRUE.
#' @param tilt.colnames Should colnames be tilted? default= FALSE.
#' @param show.row.dendrogram Plot row dendrogram? Default= TRUE.
#' @param show.row.clusters Plot row cluster names/lines? By default, they will only be showed if different clusters exist.
#' @param row.clusters.pos One of "left" or "right".
#' @param row.clusters.colors Row clusters colors (default is to use grey palette).
#' @param row.cluster.line.color Color of the line separating row clusters.
#' @param show.col.dendrogram Plot columns dendrogram? Default= TRUE.
#' @param show.col.clusters Plot column cluster names/lines?  By default, they will only be showed if different clusters exist.
#' @param col.clusters.colors Column cluster colors (default is to use grey palette).
#' @param show.legend Should the legend be plotted?
#' @param legend.title Title of the legend. Default= NA (no title).
#' @param legend.cex cex expansion factor applying to the whole legend.
#' @param display.numbers Display numbers on heatmap? Default= FALSE.
#' @param display.numbers.matrix Matrix of numbers to be displayed. By default, the x matrix will be used.
#' @param display.numbers.cex cex display numbers.
#' @param display.numbers.FUN Function to be applies to display.numbers.matrix before plotting. Default= function(x) round(x, 2).
#' @param box.lwd Line width of the box around the heatmap.
#' 
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
#' # Cluster using hierarhical clustering
#' par(mfrow=c(2,2))
#' vl_heatmap(test,
#'            cutree.rows = 3,
#'            cutree.cols = 2,
#'            show.col.clusters = F,
#'            legend.title = "log2FoldChange",
#'            main= "test")

#' # Cluster using kmeans clustering
#' res <- vl_heatmap(test,
#'                   kmeans.k = 4,
#'                   cutree.rows = 3,
#'                   cutree.cols = 5,
#'                   legend.title = "log2FoldChange",
#'                   main= "test")
#'
#' # Play with plotting parameters
#' plot(res) # Same parameters as the plot from the original function
#' plot(res,
#'      col= c("black", "blue", "yellow"),
#'      display.numbers= T,
#'      tilt.colnames= T,
#'      col.clusters.colors= rainbow(9))
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
                              clustering.distance.rows= "euclidean",
                              cutree.rows = 1L,
                              row.clusters= 1L,
                              kmeans.k= NA,
                              cluster.cols = TRUE,
                              clustering.distance.cols = "euclidean",
                              cutree.cols = 1,
                              col.clusters= 1,
                              clustering.method = "complete",
                              plot= TRUE,
                              breaks= seq(min(x, na.rm= T), max(x, na.rm= TRUE), length.out= length(col)),
                              col= c("cornflowerblue", "white", "red"),
                              na.col= "lightgrey",
                              useRaster= TRUE,
                              grid= FALSE,
                              main= NA,
                              show.rownames= TRUE,
                              show.colnames= TRUE,
                              tilt.colnames= FALSE,
                              show.row.dendrogram= TRUE,
                              show.row.clusters= length(unique(rows$cl))>1,
                              row.clusters.pos= "right",
                              row.clusters.colors= grDevices::gray.colors(100),
                              row.cluster.line.color= "black",
                              show.col.dendrogram= TRUE,
                              show.col.clusters= length(unique(cols$cl))>1,
                              col.clusters.colors= grDevices::gray.colors(100),
                              show.legend= TRUE,
                              legend.title= NA,
                              legend.cex= 1,
                              display.numbers= FALSE,
                              display.numbers.matrix= x,
                              display.numbers.cex= 1,
                              display.numbers.FUN= function(x) round(x, 2), 
                              box.lwd= .25)
{
  # Checks and default values ----
  if(is.null(rownames(x)))
    rownames(x) <- seq(nrow(x))
  if(is.null(colnames(x)))
    colnames(x) <- seq(ncol(x))
  if(!is.factor(row.clusters))
    row.clusters <- factor(row.clusters)
  if(!is.factor(col.clusters))
    col.clusters <- factor(col.clusters)
  
  # Clustering rows ----
  # Initiate rows DT
  rows <- data.table(name= rownames(x),
                     cl= row.clusters)
  rows[, order:= order(cl)]
  rcl <- NULL
  rdend <- NULL
  # Clustering
  if(cluster.rows && nrow(x)>1 && length(levels(row.clusters))==1)
  {
    set.seed(3453)
    if(is.na(kmeans.k))
    {
      .d <- if(clustering.distance.rows %in% c("pearson", "spearman"))
        as.dist(1 - cor(t(x), 
                        use= "pairwise.complete.obs", 
                        method= clustering.distance.rows)) else
                          dist(x, method = clustering.distance.rows)
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
  # Add y position
  rows[(order), y:= rev(.I)]
  
  
  # Clustering cols ----
  # Initiate cols DT
  cols <- data.table(name= colnames(x),
                     cl= col.clusters)
  cols[, order:= order(cl)]
  ccl <- NULL
  cdend <- NULL
  # Clustering
  if(cluster.cols && ncol(x)>1 && length(levels(col.clusters))==1)
  {
    set.seed(3453)
    .d <- if(clustering.distance.cols %in% c("pearson", "spearman"))
      as.dist(1 - cor(x, 
                      use= "pairwise.complete.obs", 
                      method= clustering.distance.cols)) else
                        dist(t(x), method = clustering.distance.cols)
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
  # Add x position
  cols[(order), x:= .I]
  
  # PLOT ----
  obj <- mget(ls())
  setattr(obj, "class", c("vl_heatmap", "list"))
  if(plot)
    plot(obj,
         breaks= breaks,
         col= col,
         na.col= na.col,
         useRaster= useRaster,
         grid= grid,
         main= main,
         show.rownames= show.rownames,
         show.colnames= show.colnames,
         tilt.colnames= tilt.colnames,
         show.row.dendrogram= show.row.dendrogram,
         show.row.clusters= show.row.clusters,
         row.clusters.pos= row.clusters.pos,
         row.clusters.colors= row.clusters.colors,
         row.cluster.line.color= row.cluster.line.color,
         show.col.dendrogram= show.col.dendrogram,
         show.col.clusters= show.col.clusters,
         col.clusters.colors= col.clusters.colors,
         show.legend= show.legend,
         legend.title= legend.title,
         legend.cex= legend.cex,
         display.numbers= display.numbers,
         display.numbers.matrix= display.numbers.matrix,
         display.numbers.cex= display.numbers.cex,
         display.numbers.FUN= display.numbers.FUN,
         box.lwd= box.lwd)
  
  # RETURN
  invisible(obj)
}

#' @export
plot.vl_heatmap <- function(obj,
                            breaks= NULL,
                            col= NULL,
                            na.col= NULL,
                            useRaster= NULL,
                            grid= NULL,
                            main= NULL,
                            show.rownames= NULL,
                            show.colnames= NULL,
                            tilt.colnames= NULL,
                            show.row.dendrogram= NULL,
                            show.row.clusters= NULL,
                            row.clusters.pos= NULL,
                            row.clusters.colors= NULL,
                            row.cluster.line.color= NULL,
                            show.col.dendrogram= NULL,
                            show.col.clusters= NULL,
                            col.clusters.colors= NULL,
                            show.legend= NULL,
                            legend.title= NULL,
                            legend.cex= NULL,
                            display.numbers= NULL,
                            display.numbers.matrix= NULL,
                            display.numbers.cex= NULL,
                            display.numbers.FUN= NULL, 
                            box.lwd= NULL)
{
  # Retrieve environment and update plotting variables----
  args <- intersect(ls(), names(obj))
  args <- args[sapply(args, function(x) !is.null(get(x)))]
  for(arg in args)
  {
    if(!identical(get(arg), obj[[arg]]))
      obj[[arg]] <- get(arg)
  }
  list2env(obj, environment())
  
  # Checks ----
  if(row.clusters.pos=="left")
    show.rownames <- F
  
  # Assign cluster colors ----
  rows[, col:= colorRampPalette(row.clusters.colors)(.NGRP)[.GRP], cl]
  cols[, col:= colorRampPalette(col.clusters.colors)(.NGRP)[.GRP], cl]
  
  # Heatmap image ----
  # Order matrix based on clustering
  x <- x[(rows$order), , drop=F]
  x <- x[, (cols$order), drop=F]
  # Color palette
  Cc <- circlize::colorRamp2(breaks, 
                             colors= col)
  # Image
  if(useRaster)
  {
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
  }else
  {
    rot <- t(apply(x, 2, rev))
    rot[rot<min(breaks)] <- min(breaks)
    rot[rot>max(breaks)] <- max(breaks)
    bg <- apply(rot, 2, function(x) ifelse(is.na(x), 0, NA))
    image(x= seq(nrow(bg)),
          y= seq(ncol(bg)),
          z= bg,
          breaks= c(-1, 1),
          col= na.col,
          axes= F,
          xlab= NA,
          ylab= NA)
    .br <- seq(min(breaks),
               max(breaks),
               length.out= 101)
    image(seq(nrow(rot)),
          seq(ncol(rot)),
          rot,
          axes= F,
          xlab= NA,
          ylab= NA, 
          col= Cc(.br[-101]+diff(.br)/2),
          breaks= .br,
          add= T)
  }
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
        stop("display.numbers.matrix should have the same dimensions as the heatmap!")
    # Apply function
    display.numbers.matrix <- display.numbers.FUN(display.numbers.matrix)
    
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
           col= row.cluster.line.color)
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
      text(par("usr")[1]-diff(grconvertX(c(0, par("mgp")[2]+.5), "line", "user")),
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
                    y = par("usr")[3]-diff(grconvertY(c(0, par("mgp")[2]+0.5), "line", "user")),
                    labels = colnames(x),
                    offset= 0,
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