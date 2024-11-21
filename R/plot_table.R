#' Table plot
#'
#' Plots a data.table as text
#'
#' @param DT data.table
#' @param cex Expansion factor for the table
#' @param wrap.text Number of characters after which text can be wrapped
#' @param header.col Color for header row
#' @param row.cols Colors for non-header rows
#' @return Plots table
#' @export
vl_plot_table <- function(DT,
                          cex= 1,
                          wrap.text= 20,
                          header.col= "grey50",
                          row.cols= c("grey70", "grey90"))
{
  frame <- as.data.table(as.list(names(DT)))
  setnames(frame, names(DT))
  frame <- rbind(frame, DT, fill= T)
  frame <- cbind(row= seq(nrow(frame)), frame)
  frame <- melt(frame, "row")
  frame[, value:= lapply(value, function(x) paste0(strwrap(x, width = wrap.text), collapse= "\n"))]
  frame[, value:= gsub("\n$", "", value)]
  
  plot.new()
  frame[, width:= max(strwidth(value, cex= cex))+strwidth("M", cex= cex)*2, variable]
  frame[, height:= max(strheight(value, cex= cex)+strheight("M", cex= cex)), row]
  xleft <- unique(frame[, .(variable, width)])
  xleft[, xleft:=  0.5-(sum(width)/2)+cumsum(data.table::shift(width, fill = 0))]
  frame[xleft, xleft:= xleft, on= "variable"]
  top <- unique(frame[, .(row, height)])
  top[, top:= 0.5+(sum(height)/2)-cumsum(data.table::shift(height, fill = 0))]
  frame[top, top:= top, on= "row"]
  frame[row>1, col:= rep_len(row.cols, .N), xleft]
  frame[row==1, col:= header.col]
  
  frame[, rect(xleft, top-height, xleft+width, top, xpd= T, col= col, border= "white")]
  frame[, text(xleft+width/2, top-height/2, value, xpd= T, cex= cex)]
  invisible(frame)
}