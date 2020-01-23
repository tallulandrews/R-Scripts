blank_plot <- function() {
	tmp <-  par("mar")
	par(mar=c(0,0,0,0))
	plot(1,1, col="white", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="", xlab="", ylab="", bty="n")
	return(tmp);
}
