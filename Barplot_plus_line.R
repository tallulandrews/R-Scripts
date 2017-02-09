barplot_w_line <- function(bars, line, bars.col="grey85", line.col="black", line.lwd=3, line.lty=1,...){
	if (!identical(dim(bars), dim(line))) {stop("bars and line must have the same dimensions")}
	locations = barplot(bars, beside=T, col = bars.col, ...)
	dimensions = dim(locations)
	line_y = c(0,rep(as.vector(rbind(line,c(0,0,0,0,0,0,0,0))), each=2))
	line_x = c(rep(1:((dimensions[1]+1)*dimensions[2]), each=2),(dimensions[1]+1)*dimensions[2]+1)
	lines(line_x,line_y, col=line.col, lwd=line.lwd, lty=line.lty)
}
