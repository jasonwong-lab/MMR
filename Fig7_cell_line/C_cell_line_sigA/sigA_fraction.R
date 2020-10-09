dev.new()
pdf("sigA_fraction.pdf",4,4)
par(mar=c(2,4,1,1),mgp=c(2,.5,0),bty="l")

#=============================
#       signature barplot
#=============================
mydata<-read.table("data/sigfit_fitting.txt",header = T,row.names = 1)
mydata$type<-c(rep("DLD-1",4),rep("HAP1",7))
cl<-c("#66c2a5","#fc8d62","#8da0cb")
ylim<-c(0,1)
mydata$type<-factor(mydata$type,c("HAP1","DLD-1"))
boxplot(SigA~type,data=mydata,main="Cell lines",ylab="SigA contribution",border="#cbd5e8",notch=F,ylim=ylim)
stripchart(SigA~type,data=mydata, method = "jitter",pch = 20, add = TRUE, col =cl , cex = 1,vertical=T)
dev.off()