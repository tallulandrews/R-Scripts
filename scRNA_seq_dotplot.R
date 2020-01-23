row_mean_aggregate <- function(mat, groups) {
	x <- split(seq(ncol(mat)), groups);
	results <- sapply(x, function(a) Matrix::rowMeans(mat[,a], na.rm=T))
	return(results);
}

expression_dotplot <- function(mat, clusters, genes=rownames(mat), size=0.25) {
    mat_sub <- mat[rownames(mat) %in% genes,]
    mat_detect <- row_mean_aggregate(mat_sub>0, as.factor(clusters));
    mat_sub[mat_sub==0] <- NA
    mat_mean <- row_mean_aggregate(mat_sub, clusters);
    mat_mean <- t(mat_mean)
    mat_detect <- t(mat_detect) 
    mat_mean[mat_detect==0] <- 0;
    exclude <- colSums(mat_detect)==0;
    mat_detect <- mat_detect[,!exclude]
    mat_mean <- mat_mean[,!exclude]

    splits_col <- seq(from=min(mat_mean)-0.000001, to=max(mat_mean), length=100)
    binned_col <- as.vector(t(apply(data.frame(mat_mean), 1, function(x){
                         as.numeric(cut(x, breaks=splits_col))}
			)))
    sizes <- as.vector(mat_detect)*size
# from scClustViz
    layout(matrix(c(0,2,3,1),2),widths=c(1,5),heights=c(1,5))
    par(mar=c(9,0,0,.5))
    plot(x=NULL,y=NULL,xlim=c(0.5,nrow(mat_mean)+.5),ylim=c(0.5,ncol(mat_mean)+.5),
         xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab=NA,ylab=NA,bty="n")
    abline(v=1:nrow(mat),col="grey90")
    symbols(x=rep(1:nrow(mat_mean),ncol(mat_mean)),
            y=as.vector(sapply(1:ncol(mat_mean),function(X) rep(X,nrow(mat_mean)))),
            circles=sizes,inches=F,add=T,xpd=NA,
            fg=NA, #colorspace::sequential_hcl(100,palette="Viridis",rev=T)[binned_col],
            bg=colorspace::sequential_hcl(100,palette="Viridis",rev=T)[binned_col])
    axis(side=2,at=1:ncol(mat_mean),lwd=0,labels=colnames(mat_mean),las=2,cex.axis=0.8)
    axis(side=1,at=1:nrow(mat_mean),lwd=0,labels=rownames(mat_mean),las=2,cex.axis=1.2)
    # Legend:
    tx0 <- par("usr")[1]
    tx <- (par("usr")[2] - par("usr")[1])
    ty0 <- par("usr")[3]
    ty <- par("usr")[4] - par("usr")[3]
    segments(x0=tx0 - seq(.15,.03,length.out=1000) * tx,
             y0=ty0 - 0.02 * ty,y1=ty0 - 0.05 * ty,
             col=colorspace::sequential_hcl(1000,palette="Viridis",rev=T),xpd=NA)
    text(x=tx0 - c(.15,.09,.03) * tx,
         y=ty0 - c(0.035,0.02,0.035) * ty,
         labels=c(round(min(splits_col),2),
                  "Avg detected expr",
                  round(max(splits_col),2)),pos=2:4,xpd=NA)
    symbols(x=tx0 - c(.15,.09,.03) * tx,
            y=ty0 - rep(.14,3) * ty,add=T,xpd=NA,
            circles=c(0.25,0.5,0.75)*size,inches=F,bg="black", fg=NA)
    text(x=tx0 - c(.149,.089,0.029,.09) * tx,
         y=ty0 - c(rep(.23,3),.26) * ty,xpd=NA,
         labels=c("25%","50%","75%","Detection Rate"))
}
