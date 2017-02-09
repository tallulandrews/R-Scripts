Density_Coloured_Scatterplot <- function (x, y, gradient_colours=c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100")) {
	
        fancy <- densCols(x, y, colramp=colorRampPalette(c("black","white")))
        dens <- col2rgb(fancy)[1,]+1L
        colours <-  colorRampPalette(gradient_colours)(256) #blue->yellow
        dens.col <- colours[dens]

     	plot(x,y, main="", ylab="", xlab="", col = dens.col,pch=16)
	invisible(dens.col);
}
