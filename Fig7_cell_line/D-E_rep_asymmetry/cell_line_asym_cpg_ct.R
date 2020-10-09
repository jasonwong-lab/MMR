dev.new()
pdf("cell_line_asym_cpg_ct.pdf",8,4)
layout(matrix(1:2,nrow=1,byrow = T))
par(mar=c(2,4,1,1),mgp=c(2,.5,0),bty="l")
ylim<-c(-1,0.5)

#=================
#   HAP1
#================
cpg<-read.table("data/HAP1/ct_cpg.txt",header = T,row.names = 1)
cpg$ct<-log(cpg$ct,2)
cpg<-as.data.frame(cpg[-8,])
colnames(cpg)<-"value"
cpg$type<-rep("cpg",nrow(cpg))
ncpg<-read.table("data/HAP1/ct_ncpg.txt",header = T,row.names = 1)
ncpg$ct<-log(ncpg$ct,2)
ncpg<-as.data.frame(ncpg[-8,])
colnames(ncpg)<-"value"
ncpg$type<-rep("ncpg",nrow(ncpg))
muts.ct1<-rbind(cpg,ncpg)
boxplot(value~type,data=muts.ct1,main="HAP1 mutants",ylab="Log2(Lagging/Leading)",border="#cbd5e8",notch=F,names=c("CpG","Non-CpG"),ylim=ylim,cex=.8)
stripchart(value~type,data=muts.ct1, method = "jitter",pch = 20, add = TRUE, col =c("#984ea3","#ff7f00") , cex = 1,vertical=T)
t<-t.test(muts.ct1$value[which(muts.ct1$type=="cpg")],muts.ct1$value[which(muts.ct1$type=="ncpg")])
mtext(paste0("P=",round(t$p.value,5)),3,cex=.7,line=-1)
abline(h=0,col="#636363",lty=2)

#=================
#   DLD1
#================
cpg<-read.table("data/DLD1/ct_cpg.txt",header = T,row.names = 1)
cpg$ct<-log(cpg$ct,2)
cpg<-as.data.frame(cpg)
colnames(cpg)<-"value"
cpg$type<-rep("cpg",nrow(cpg))
ncpg<-read.table("data/DLD1/ct_ncpg.txt",header = T,row.names = 1)
ncpg$ct<-log(ncpg$ct,2)
ncpg<-as.data.frame(ncpg)
colnames(ncpg)<-"value"
ncpg$type<-rep("ncpg",nrow(ncpg))
muts.ct<-rbind(cpg,ncpg)
boxplot(value~type,data=muts.ct,main="DLD-1 mutants",ylab="Log2(Lagging/Leading)",border="#cbd5e8",notch=F,names=c("CpG","Non-CpG"),ylim=ylim)
stripchart(value~type,data=muts.ct, method = "jitter",pch = 20, add = TRUE, col =c("#984ea3","#ff7f00") , cex = 1,vertical=T)
t<-t.test(muts.ct$value[which(muts.ct$type=="cpg")],muts.ct$value[which(muts.ct$type=="ncpg")])
mtext(paste0("P=",signif(t$p.value,3)),3,cex=.7,line=-1)
abline(h=0,col="#636363",lty=2)

dev.off()