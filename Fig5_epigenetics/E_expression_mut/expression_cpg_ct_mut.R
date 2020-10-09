dev.new()
pdf("expression_cpg_ct_mut.pdf",4,4)
#layout(matrix(1:4,nrow=1))
par(mar=c(2,4,1,1),mgp=c(2,.5,0),bty="l")

#=========================
#  expression and mutation
#=========================
#============MBD4=================
obs<-read.table("data/MBD4/obs_mut.txt",row.names = 1)
exp<-read.table("data/MBD4/exp_mut.txt",row.names = 1)
ratio.mbd4.cpg<-as.numeric(obs[1,][2:5])/as.numeric(exp[1,][1:4])

#============MSI=================
obs<-read.table("data/MSI/obs_mut.txt",row.names = 1)
exp<-read.table("data/MSI/exp_mut.txt",row.names = 1)
ratio.msi.cpg<-as.numeric(obs[1,][2:5])/as.numeric(exp[1,][1:4])

#============MSS=================
obs<-read.table("data/MSS/obs_mut.txt",row.names = 1)
exp<-read.table("data/MSS/exp_mut.txt",row.names = 1)
ratio.mss.cpg<-as.numeric(obs[1,][2:5])/as.numeric(exp[1,][1:4])

#-------------------------plot--------------------------
x<-1:4
cl<-c("#ef8a62","#67a9cf","#bdbdbd")

#------------cpg c>t slope-------
y<-ratio.mbd4.cpg
y2<-ratio.mss.cpg
y3<-ratio.msi.cpg

ylim<-range(y,y2,y3)
par(mar=c(3.2,4,0.8,1),mgp=c(2,.5,0))
plot(x,y,xlab="Expression(Low->High)",ylab="Ratio(Obs/Exp)",bty="l",ylim=ylim,pch=19,col=cl[1],main="CpG C>T")
lines(as.numeric(y),lty=2,col=cl[1])
points(x,y2,pch=19,lty=2,col=cl[2])
lines(as.numeric(y2),lty=2,col=cl[2])

points(x,y3,pch=19,lty=2,col=cl[3])
lines(as.numeric(y3),lty=2,col=cl[3])

legend("topright",legend=c("MBD4 mutants","MSS","MMRd(MSI)"),pch=c(19,19,19),col=cl[1:3],xpd=T,box.col=NA,pt.cex=1,cex=1,ncol=1,x.intersp=.2,lty=2,y.intersp=.8)
dev.off()