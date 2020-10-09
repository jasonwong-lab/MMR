dev.new()
cairo_pdf("h3k36me3_mut.pdf",3,3)
par(mar=c(4,4,1,1),mgp=c(2,.5,0))

#==========================
#    plot h3 and mutation
#==========================
mydata<-read.table("data/methyl_bin_mut.txt")
h3<-aggregate(V3~V1, data = mydata, mean)
rownames(h3)<-h3$V1
bins<-paste0("bin",1:10)
h3<-h3[bins,]
x<-1:10
y1<-h3$V3
cl<-brewer.pal(4,"Set1")
l.col<-c()
t<-80
l.col[1]<-rgb(228,26,28,t,maxColorValue=228)
l.col[2]<-rgb(55,126,184,t,maxColorValue=184)

plot(x,y1,xlab="H3K36me3 bins (Low->High)",ylab="Ratio(obs/exp)",bty="l",col=cl[1],pch=19,lty=2)
fit<-lm(y1~x)
abline(fit,lty=2,col=l.col[1])
tmp<-cor.test(x,y1)
r<-signif(tmp$estimate,3)
alpha<-signif(fit$coefficients[2],3)
mtext(paste0("\u03b1","=",alpha,", P<0.0001"),3,cex=.7,line=-1)
dev.off()
