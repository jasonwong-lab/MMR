pdf("venn_diagram.pdf",4,4)
cl<-c("#66c2a5","#fc8d62","#8da0cb","#e78ac3")
a<-c(141,160,203)
b<-c(231,138,195)
cl.ol<-rgb(372,298,398,maxColorValue=398)
par(mar=c(1,1,1,1))
plot(0,type="n",bty="n",ann=F,axes=F,xlim=c(-2.5,2.5),ylim=c(-2.5,2.5))
r=1.2
draw.circle(-1,1,r,col=cl[1])
draw.circle(1,1,r,col=cl[2])
draw.circle(-1,-1,r,col=cl[3])
draw.circle(1,-1,r,col=cl[4])

# draw overlap
# top
y1<-sqrt(r^2-1)+1
y2<-1-sqrt(r^2-1)
yvals <- seq(y2, y1, 0.01)
xvals <- sqrt(r^2 - (yvals-1)^2) - 1
yvals <- c(yvals, rev(yvals))
xvals <- c(xvals, -xvals)
polygon(xvals,yvals,col="#f5f5f5")

# bot
y1<-sqrt(r^2-1)+1
y2<-1-sqrt(r^2-1)
yvals <- seq(y2, y1, 0.01)
xvals <- sqrt(r^2 - (yvals-1)^2) - 1
yvals <- c(-yvals, -rev(yvals))
xvals <- c(xvals, -xvals)
polygon(xvals,yvals,col=cl.ol)

#left
y1<-sqrt(r^2-1)+1
y2<-1-sqrt(r^2-1)
yvals <- seq(y2, y1, 0.01)
xvals <- sqrt(r^2 - (yvals-1)^2) - 1
yvals <- c(-yvals, -rev(yvals))
xvals <- c(xvals, -xvals)
polygon(yvals,xvals,col=cl[1])

#right
y1<-sqrt(r^2-1)+1
y2<-1-sqrt(r^2-1)
yvals <- seq(y2, y1, 0.01)
xvals <- sqrt(r^2 - (yvals-1)^2) - 1
yvals <- c(yvals, rev(yvals))
xvals <- c(xvals, -xvals)
polygon(yvals,xvals,col=cl[2])

dev.off()