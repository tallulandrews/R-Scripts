color.bar <- function(lut, min=0, max=1, ticks.at=seq(min, max,len=min(10,length(lut)+1)), ticks.lab=ticks_at, title='', horiz=FALSE) {
    scale = (length(lut))/(max-min)

    if (horiz) {
        dev.new(width=5, height=1.75)
        par(mar=c(4,0,1,0))
        plot(c(min,max),c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
        axis(1, at=ticks.at, labels=ticks.lab, las=1)
        for (i in 1:(length(lut))) {
             x = (i-1)/scale + min
             rect(x,0,x+1/scale,10, col=lut[i], border=NA)
        }
        mtext(side=1, line=2.5, title, font=2, cex=1.1)

    } else {
        dev.new(width=1.75, height=5)
        plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
        axis(2, at=ticks.at, labels=ticks.lab, las=1)
        for (i in 1:(length(lut))) {
             y = (i-1)/scale + min
             rect(0,y,10,y+1/scale, col=lut[i], border=NA)
        }
    }
}

