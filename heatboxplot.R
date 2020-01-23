silhouette_boxplot <- function(x, col="grey50", ylim=c(min(x),max(x)), main="") {
	par(mar=c(0,0,1,0))
	boxplot(x, ylim=ylim, col=col, border=col, staplewex=0, outline=FALSE, xaxt="n", yaxt="n", xlab="", ylab="", frame.plot=FALSE)

}
colorbar_plot <- function(colour_vec, horiz=FALSE) {
	par(mar=c(1,0,1,0))
	if (!horiz) {
		image(rbind(1:length(colour_vec)), 
			col = colour_vec, axes = FALSE)
	} else {
		image(cbind(1:length(colour_vec)), 
			col = colour_vec, axes = FALSE)
	}
}

get_height_for_k <- function(dendro, k) {
	height_out <- -1;
        curr_k <- 1;
        dendro_list <- list(dendro)
        dendro_heights <- attr(dendro, "height")
        while( curr_k < k ){
                to_split <- which(dendro_heights == max(dendro_heights))
                to_split_dendro <- dendro_list[[to_split]]
                to_split_height <-  dendro_heights[to_split]

                children <- as.list(to_split_dendro)
                for (i in 1:length(children)) {
                        dendro_heights <- c(dendro_heights,attr(children[[i]],"height"))
                        dendro_list[[length(dendro_list)+1]] <- children[[i]]
                }
                # Remove to split
                dendro_list[to_split] <- NULL
		height_out = mean(c(dendro_heights[to_split], 
				max(dendro_heights[-to_split])));
                dendro_heights <- dendro_heights[-to_split]
                curr_k <- curr_k-1+length(children)
        }
	return(height_out)
}



cut_dendro_plot <- function(dendro, k, horiz=TRUE, full=TRUE) {
	par(mar=c(0,0,0,0))
	dendro <- as.dendrogram(dendro);
	if (!full) { dendro <- cut(dendro, k=k) }
	plot(dendro, axes=FALSE, leaflab="none", horiz=horiz)

	if (full) {
		h <- get_height_for_k(dendro, k)
		abline(v=h, lty=2, col="grey40")
	}
}



epic_dendro_boxplots <- function(x, distfun = dist, hclustfun = hclust,
	k, markers, marker_col, vbar_labels, vbar_cols, fixed_ylim=NA) {
	# Layout
	box_height=1; box_width=0.5;
	bar_width=0.25; dendro_width=3;

	nrow = k; ncol = length(markers);
        lhei <- c(rep(box_height, times=k))
        lwid <- c(dendro_width, bar_width, rep(box_width, times=ncol), axis_width)
        lmat <- matrix(c(rep(1,times=k),rep(2,times=k),3:(ncol*k+2), rep(ncol*k+3, times=k)), 
			ncol=(ncol+3), byrow=FALSE)
print(dim(lmat))

	layout(lmat, heights=lhei, widths=lwid, respect=FALSE)

	# Create Dendrogram
        hcc <- hclustfun(distfun( t(x) ))
	cut_dendro_plot(hcc, k)

	# Create colour bar
        ddc <- as.dendrogram(hcc)
#       ddc <- reorder(ddc, colnames(x))
        colInd <- order.dendrogram(ddc)
	if (!is.factor(vbar_labels)) {
		vbar_labels <- factor(vbar_labels);
	}
	my_col <- vbar_cols[vbar_labels[colInd]]
	colorbar_plot(my_col)

	# Create marker boxes
	groups <- cutree(as.hclust(ddc), k=k)
	groups <- groups[colInd]
	x <- x[,colInd]
	g_order <- rev(unique(groups))
	for (m in markers) {
		if (sum(is.na(fixed_ylim))) {
			ylim = c(min(x[rownames(x) == m,]), max(x[rownames(x) == m,]))
		} else {
			ylim = fixed_ylim
		}
		global = mean(x[rownames(x) == m,]);
		for (g in g_order) {
			dat <- x[rownames(x) == m, groups == g]
			silhouette_boxplot(dat, col=marker_col[markers==m], ylim=ylim)
			points(1,global, col="black", pch=18, cex=1.5)
			if (g == g_order[1]) {
				title(main=m)
			}
			if (m == markers[length(markers)] & !sum(is.na(fixed_ylim))) {
				par(xpd=TRUE)
				axis(4, line=-1)
			}
		}
	}
}
