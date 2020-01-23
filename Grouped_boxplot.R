grouped_boxplot <- function(x, y, z, cols=brewer.pal(length(unique(z)), "Greys"), place_legend=NULL, legend_title=NULL, ...) {
        nx <- length(levels(x))
        ny <- length(levels(y))
        nz <- length(levels(z))
        xes <- 1:(nx*(nz+1))
        xes <- xes[-seq(from=nz+1, to=length(xes), by=nz+1)]
        thing <- boxplot(y~z+x, at=xes, col=cols, xaxt="n", names=rep("", length(xes)), ...)
        x_names <- rep(levels(x), each=nz)
        x_axis_info <- aggregate(xes, by=list(x_names), mean)
        axis(1, at=x_axis_info[,2], labels=x_axis_info[,1], las=2)
        if (!is.null(place_legend)) {
                legend(place_legend, bty="n", fill=cols, levels(z), title=legend_title)
        }
	invisible(thing)
}
