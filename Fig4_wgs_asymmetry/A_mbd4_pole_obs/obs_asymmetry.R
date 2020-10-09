muts<-"MBD4" #n=18
mutl<-"TCGA-F5-6814" #n=197

p3<-read.table("data/prime3.txt",row.names = 1)
p5<-read.table("data/prime5.txt",row.names = 1)

#===============get muts data==================
muts.p5<-apply(p5[muts,],2,sum)
muts.p3<-apply(p3[muts,],2,sum)

#===============get mutl data==================
mutl.p5<-apply(p5[mutl,],2,sum)
mutl.p3<-apply(p3[mutl,],2,sum)

muts.data1<-c(muts.p5,muts.p3)
mutl.data1<-c(mutl.p5,mutl.p3)

order<-c(7,8,9,10,11,12,6,5,4,3,2,1,19,20,21,22,23,24,18,17,16,15,14,13)
muts.data<-muts.data1[order]
mutl.data<-mutl.data1[order]

dev.new()
pdf("obs_rep_asym.pdf",4,4)
layout(matrix(1:3,nrow=3),heights = c(1,4,4))
cl1<-c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
cl<-rep(cl1,each=2)

###plot legend
par(mar=c(0,2,1,2),mgp=c(1.7,.5,.5))
plot(0,type="n",bty="n",ann=F,axes=F)
legend("center",legend=c("C>A/G>T","C>G/G>C","C>T/G>A","T>A/A>T","T>C/A>G","T>G/A>C"),col=cl1,border=cl1,xpd=T,box.col=NA,pch=15,pt.cex=2,cex=.8,ncol=6,text.width=.1,y.intersp=2)

par(mar=c(2,4,1,1),mgp=c(2,.5,.5))
barplot(muts.data,beside=F,col=cl,border=cl,space= rep(rep(c(1,0.2),6),2),axisnames=F,axes=T,ylab="")
mtext("MBD4 mutants",side=3,adj=0)
barplot(mutl.data,beside=F,col=cl,border=cl,space= rep(rep(c(1,0.2),6),2),axisnames=F,axes=T,ylab="")
mtext("POLE mutants",side=3,adj=0)
dev.off()
