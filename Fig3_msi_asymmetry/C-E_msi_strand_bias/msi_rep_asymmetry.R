dev.new()
cairo_pdf("msi_rep_asymmetry.pdf",8,2)
layout(matrix(c(1,2,3),nrow=1),widths = c(1,1,1.5))
par(mar=c(4,4,1,1),mgp=c(2,.5,0),bty="l")
cl<-c("#67a9cf","#ef8a62")
ylim=c(-2,1)

#=================
#   PLOT T>C
#=================
p3_all<-read.table("data/mut_prime3.txt",row.names = 1)
p5_all<-read.table("data/mut_prime5.txt",row.names = 1)
share_name<-intersect(rownames(p3_all),rownames(p5_all))
p3_all<-p3_all[share_name,]
p5_all<-p5_all[share_name,]
rep.all<-cbind(p3_all,p5_all)
#-----------------read muts-----------
msi.all<-read.table("data/msi_anno.txt",header = T,row.names = 1)
msi.all<-subset(msi.all,muts!="Double")
msi.all$muts<-as.character(msi.all$muts)
#---------tc rm zero--------------
rep.all$agall<-rep.all[,3]+rep.all[,16]
rep.all$tcall<-rep.all[,4]+rep.all[,15]
p3<-subset(rep.all,select=c("agall","tcall"))
for (i in 1: nrow(p3)){
  p3[i,][which(p3[i,]<4)]=NA
}
p3.rm<-na.omit(p3)
n<-nrow(p3.rm)
p3.rm.ratio<-as.data.frame(matrix(rep(NA,n),ncol=1))
p3.rm.ratio[,1]<-p3.rm[,1]/p3.rm[,2]
rownames(p3.rm.ratio)<-rownames(p3.rm)
colnames(p3.rm.ratio)<-"tc"
msi.all1<-msi.all[rownames(p3.rm),]
msi.all1<-cbind(p3.rm.ratio,msi.all1)
mmr.all1<-msi.all1
mmr.all1$tc<-log(mmr.all1$tc,2)
tmp<-data.frame(tc=mmr.all1$tc)
rownames(tmp)<-rownames(mmr.all1)
boxplot(tc~muts,data=mmr.all1,ylab="Log2(Lagging/Leading)",border="#cbd5e8",notch=F,names=c(paste0("MutL","\u03b1"),paste0("MutS","\u03b1")),main="A>G/T>C",ylim=ylim)
stripchart(tc~muts,data=mmr.all1, method = "jitter",pch = 20, add = TRUE, col =cl , cex = 1,vertical=T)
abline(h=0,lty=2,col="#bdbdbd")
t<-t.test(mmr.all1$tc[which(mmr.all1$muts==1)],mmr.all1$tc[which(mmr.all1$muts==0)])
mtext(paste0("P=",signif(t$p.value,3)),3,cex=.7,line=-1)


#=================
#   PLOT C>T
#=================
#---------ct rm zero--------------
rep.all$ctall<-rep.all[,11]+rep.all[,24]
rep.all$gaall<-rep.all[,12]+rep.all[,23]
p3<-subset(rep.all,select=c("ctall","gaall"))
for (i in 1: nrow(p3)){
  p3[i,][which(p3[i,]<4)]=NA
}
p3.rm<-na.omit(p3)
n<-nrow(p3.rm)
p3.rm.ratio<-as.data.frame(matrix(rep(NA,n),ncol=1))
p3.rm.ratio[,1]<-p3.rm[,1]/p3.rm[,2]
rownames(p3.rm.ratio)<-rownames(p3.rm)
colnames(p3.rm.ratio)<-"ct"
msi.all1<-msi.all[rownames(p3.rm),]
msi.all1<-cbind(p3.rm.ratio,msi.all1)
mmr.all1<-msi.all1
mmr.all1$ct<-log(mmr.all1$ct,2)
boxplot(ct~muts,data=mmr.all1,ylab="Log2(Lagging/Leading)",border="#cbd5e8",notch=F,names=c(paste0("MutL","\u03b1"),paste0("MutS","\u03b1")),main="C>T/G>A",ylim=ylim)
stripchart(ct~muts,data=mmr.all1, method = "jitter",pch = 20, add = TRUE, col =cl , cex = 1,vertical=T)
abline(h=0,lty=2,col="#bdbdbd")
t<-t.test(mmr.all1$ct[which(mmr.all1$muts==1)],mmr.all1$ct[which(mmr.all1$muts==0)])
mtext(paste0("P=",signif(t$p.value,3)),3,cex=.7,line=-1)


#=====================================================================
#   MutS mutants
#======================================================================
ylim<-c(-1.6,1.5)
cpg<-read.table("data/muts_ct_cpg.txt",header = T,row.names = 1)
colnames(cpg)<-"value"
cpg$type<-rep("cpg",nrow(cpg))
ncpg<-read.table("data/muts_ct_ncpg.txt",header = T,row.names = 1)
colnames(ncpg)<-"value"
ncpg$type<-rep("ncpg",nrow(ncpg))
muts.ct1<-rbind(cpg,ncpg)


#=================
#   MutL mutants
#================
cpg<-read.table("data/mutl_ct_cpg.txt",header = T,row.names = 1)
colnames(cpg)<-"value"
cpg$type<-rep("cpg1",nrow(cpg))
ncpg<-read.table("data/mutl_ct_ncpg.txt",header = T,row.names = 1)
colnames(ncpg)<-"value"
ncpg$type<-rep("ncpg1",nrow(ncpg))
muts.ct<-rbind(cpg,ncpg)
muts.ct.msi<-rbind(muts.ct1,muts.ct)
muts.ct.msi$type<-factor(muts.ct.msi$type,c("cpg1","cpg","ncpg1","ncpg"))

par(mar=c(4,4,1,1),mgp=c(2,.5,0),bty="l")
boxplot(value~type,data=muts.ct.msi,main="MMRd (WXS)",ylab="Log2(Lagging/Leading)",border="#cbd5e8",notch=F,names=c(paste0("MutL","\u03b1"),paste0("MutS","\u03b1"),paste0("MutL","\u03b1"),paste0("MutS","\u03b1")),ylim=ylim,cex.axis=1,width=c(0.5,0.5,0.5,0.5))
stripchart(value~type,data=muts.ct.msi, method = "jitter",pch = 20, add = TRUE, col =cl, cex = 1,vertical=T)
abline(h=0,lty=2,col="#bdbdbd")
mtext("CpG",1,cex=.8,line=2,adj=.25)
mtext("Non-CpG",1,cex=.8,line=2,adj=.85)
dev.off()

