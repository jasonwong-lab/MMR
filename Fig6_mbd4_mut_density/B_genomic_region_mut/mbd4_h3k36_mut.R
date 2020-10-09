#-----------------------------------
#      obs/exp
#-----------------------------------
a<-0
mbd4.high<-read.table("data/top20MBD4_notop20H3K36_met.bed")
mbd4.high<-na.omit(mbd4.high)
mbd4.high<-subset(mbd4.high,V26>a)

h3k36.high<-read.table("data/top20H3K36_notop20MBD4_met.bed")
h3k36.high<-na.omit(h3k36.high)
h3k36.high<-subset(h3k36.high,V26>a)

mbd4.low<-read.table("data/bot20MBD4_notop20H3K36_met.bed")
mbd4.low<-na.omit(mbd4.low)
mbd4.low<-subset(mbd4.low,V26>a)

h3k36.low<-read.table("data/bot20H3K36_notop20MBD4_met.bed")
h3k36.low<-na.omit(h3k36.low)
h3k36.low<-subset(h3k36.low,V26>a)
data.list<-list(mbd4.high,h3k36.high,mbd4.low,h3k36.low)
obs.exp<-data.frame(mbd4=rep(NA,4),mss=rep(NA,4),msi=rep(NA,4))
for(i in 1:4){
  tmp<-data.list[[i]]
  tmp1<-sum(tmp$V8)
  tmp2<-sum(tmp$V9)
  obs.exp[i,1]<-tmp1/tmp2
  
  tmp1<-sum(tmp$V10)
  tmp2<-sum(tmp$V11)
  obs.exp[i,2]<-tmp1/tmp2
  
  tmp1<-sum(tmp$V12)
  tmp2<-sum(tmp$V13)
  obs.exp[i,3]<-tmp1/tmp2
}

#  plot
dev.new()
pdf("mbd4_h3k36_mut.pdf",4,4)
par(mar=c(3,3,1,1),mgp=c(1.5,.5,0),bty="o")

cl<-c("#ffeda0","#feb24c","#f03b20")
at<-barplot(t(obs.exp),beside=T,col=cl,axes=T,ylab="Obs/Exp",axisnames=F,ylim = c(0,1.5))
axis(1,at=(at[1,]+at[2,]+at[3,])/3,labels = c("MBD4 high","H3K36me3 high","MBD4 low","H3K36me3 low"),cex.axis=.7)
abline(h=1,col="#636363",lty=2)
legend(0.4,1.5,legend=c("MBD4 mutants","MSS","MMRd(MSI)"),fill=cl,xpd=T, box.col=NA,ncol=1,x.intersp = .5)
dev.off()
