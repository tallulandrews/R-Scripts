rotate<-function(x) t(apply(x,2,rev))
fancy_table_plot <- function(mat) {
        names=colnames(mat)
        greens = c("#e5f5e0","#a1d99b","#31a354")
        blues = c("#deebf7","#9ecae1","#2182bd")
        purples = c("#efedf5","#bcbddc","#756bb1")

        my_palette = c(rev(greens),"white","#fc9272","grey85","#fc9272","white",blues,purples)
        mat_col=mat
        for (i in 1:length(names)) {
        for (j in 1:length(names)) {
                if (i > j) {
                        mat_col[i,j] = -mat_col[i,j]
                }
                if (i == j) {
                        if(mat_col[i,j] != 0) { mat_col[i,j] = mat_col[i,j]+100}
                }
        }
        }
        my_breaks = c(-100,-2,-1.2,-1,-0.8,-0.000001, 0.0000001, 0.8,1,1.2,2,100,101.2,102,200)

        par(mar=c(1,4.5,4.5,1))
        image(rotate(mat_col), col=my_palette, breaks=my_breaks, xaxt="n",yaxt="n")
        box()
        coords = seq(from=0, to=1, length=length(names))
        axis(2,at=coords,labels=rev(names),tick=FALSE,las=2)
        axis(3,at=coords, labels=names, tick=FALSE)
        cell_contents = round(as.vector(mat),digits=3)
        cell_contents = as.character(cell_contents)
        cell_contents[cell_contents == "0"] = "NA"
        text(rep(coords,each=length(names)), rep(rev(coords),times=length(names)), cell_contents)
}

