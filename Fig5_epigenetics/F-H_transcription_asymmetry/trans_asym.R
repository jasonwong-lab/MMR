dev.new()
pdf("trans_asym.pdf",6,2)
layout(matrix(1:3,nrow=1))
par(mar=c(2,4,1,1),mgp=c(2,.5,0),bty="l")

ylim<-c(-0.4,0.4)
#==================
#    MSI Mutants
#==================
cpg<-read.table("data/MSI/ct_cpg.txt",header = T,row.names = 1)
cpg$ct<-log(cpg$ct,2)
cpg<-as.data.frame(cpg[-8,])
colnames(cpg)<-"value"
cpg$type<-rep("cpg",nrow(cpg))
ncpg<-read.table("data/MSI/ct_ncpg.txt",header = T,row.names = 1)

ncpg$ct<-log(ncpg$ct,2)
ncpg<-as.data.frame(ncpg[-8,])
colnames(ncpg)<-"value"
ncpg$type<-rep("ncpg",nrow(ncpg))
muts.ct<-rbind(cpg,ncpg)
boxplot(value~type,data=muts.ct,main="MMRd(MSI)",ylab="Log2(Template/Coding)",border="#cbd5e8",notch=F,names=c("CpG","Non-CpG"),ylim=ylim)
stripchart(value~type,data=muts.ct, method = "jitter",pch = 20, add = TRUE, col =c("#984ea3","#ff7f00") , cex = 1,vertical=T)
t<-t.test(muts.ct$value[which(muts.ct$type=="cpg")],muts.ct$value[which(muts.ct$type=="ncpg")])
mtext(paste0("P=",signif(t$p.value,3)),3,cex=.7,line=-1)

#===================
# MBD4 mutants
#===================
cpg<-read.table("data/MBD4/ct_cpg.txt",header = T,row.names = 1)
cpg$ct<-log(cpg$ct,2)
cpg<-as.data.frame(cpg)
colnames(cpg)<-"value"
cpg$type<-rep("cpg",nrow(cpg))
ncpg<-read.table("data/MBD4/ct_ncpg.txt",header = T,row.names = 1)

ncpg$ct<-log(ncpg$ct,2)
ncpg<-as.data.frame(ncpg)
colnames(ncpg)<-"value"
ncpg$type<-rep("ncpg",nrow(ncpg))
muts.ct<-rbind(cpg,ncpg)
# muts.ct<-rbind(muts.ct,muts.ct1)

boxplot(value~type,data=muts.ct,main="MBD4 mutants",ylab="Log2(Template/Coding)",border="#cbd5e8",notch=F,names=c("CpG","Non-CpG"),ylim=ylim)
stripchart(value~type,data=muts.ct, method = "jitter",pch = 20, add = TRUE, col =c("#984ea3","#ff7f00") , cex = 1,vertical=T)
t<-t.test(muts.ct$value[which(muts.ct$type=="cpg")],muts.ct$value[which(muts.ct$type=="ncpg")])
mtext(paste0("P=",signif(t$p.value,3)),3,cex=.7,line=-1)

#==================
#    MSS Mutants
#==================
cpg<-read.table("data/MSS/ct_cpg.txt",header = T,row.names = 1)
cpg$ct<-log(cpg$ct,2)
cpg<-as.data.frame(cpg)
colnames(cpg)<-"value"
cpg$type<-rep("cpg",nrow(cpg))
ncpg<-read.table("data/MSS/ct_ncpg.txt",header = T,row.names = 1)

ncpg$ct<-log(ncpg$ct,2)
ncpg<-as.data.frame(ncpg)
colnames(ncpg)<-"value"
ncpg$type<-rep("ncpg",nrow(ncpg))

muts.ct<-rbind(cpg,ncpg)
boxplot(value~type,data=muts.ct,main="MSS",ylab="Log2(Template/Coding)",border="#cbd5e8",notch=F,names=c("CpG","Non-CpG"),ylim=ylim)
stripchart(value~type,data=muts.ct, method = "jitter",pch = 20, add = TRUE, col =c("#984ea3","#ff7f00") , cex = 1,vertical=T)
t<-t.test(muts.ct$value[which(muts.ct$type=="cpg")],muts.ct$value[which(muts.ct$type=="ncpg")])
mtext(paste0("P=",signif(t$p.value,3)),3,cex=.7,line=-1)

dev.off()