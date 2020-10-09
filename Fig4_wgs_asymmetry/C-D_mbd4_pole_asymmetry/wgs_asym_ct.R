dev.new()
pdf("wgs_asym_ct.pdf",8,4)
layout(matrix(1:2,nrow=1,byrow = T))
par(mar=c(2,3,1,1),mgp=c(2,.5,0),bty="l")
ylim<-c(-1,1)
#=================
#   MBD4
#================
cpg<-read.table("data/mbd4_ct_cpg.txt",header = T,row.names = 1)
cpg$ct<-log(cpg$ct,2)
cpg<-as.data.frame(cpg)
colnames(cpg)<-"value"
cpg$type<-rep("cpg",nrow(cpg))
ncpg<-read.table("data/mbd4_ct_ncpg.txt",header = T,row.names = 1)
ncpg$ct<-log(ncpg$ct,2)
ncpg<-as.data.frame(ncpg)

colnames(ncpg)<-"value"
ncpg$type<-rep("ncpg",nrow(ncpg))
muts.ct<-rbind(cpg,ncpg)

boxplot(value~type,data=muts.ct,main="MBD4 mutants (WGS)",ylab="Log2(Lagging/Leading)",border="#cbd5e8",notch=F,names=c("CpG","Non-CpG"),ylim=ylim)
stripchart(value~type,data=muts.ct, method = "jitter",pch = 20, add = TRUE, col =c("#984ea3","#ff7f00") , cex = 1,vertical=T)
abline(h=0,lty=2,col="#bdbdbd")
t<-t.test(muts.ct$value[which(muts.ct$type=="cpg")],muts.ct$value[which(muts.ct$type=="ncpg")])
mtext(paste0("P=",signif(t$p.value,3)),3,cex=.7,line=-1)

#=================
#   POLE
#================
sap<-c("TCGA-A6-6141","TCGA-AA-3555","TCGA-AA-3977","TCGA-AA-A00N","TCGA-AZ-4315","TCGA-CA-6717","TCGA-CA-6718","TCGA-EI-6917","TCGA-F5-6814")
cpg<-read.table("data/pole_ct_cpg.txt",header = T,row.names = 1)
cpg$ct<-log(cpg$ct,2)
cpg<-as.data.frame(cpg)
colnames(cpg)<-"value"
cpg$type<-rep("cpg",nrow(cpg))
cpg<-cpg[sap,]
ncpg<-read.table("data/pole_ct_ncpg.txt",header = T,row.names = 1)
ncpg$ct<-log(ncpg$ct,2)
ncpg<-as.data.frame(ncpg)
colnames(ncpg)<-"value"
ncpg$type<-rep("ncpg",nrow(ncpg))
ncpg<-ncpg[sap,]
muts.ct<-rbind(cpg,ncpg)
boxplot(value~type,data=muts.ct,main="POLE mutants (WGS)",ylab="Log2(Lagging/Leading)",border="#cbd5e8",notch=F,names=c("CpG","Non-CpG"),ylim=ylim)
stripchart(value~type,data=muts.ct, method = "jitter",pch = 20, add = TRUE, col =c("#984ea3","#ff7f00") , cex = 1,vertical=T)
abline(h=0,lty=2,col="#bdbdbd")
t<-t.test(muts.ct$value[which(muts.ct$type=="cpg")],muts.ct$value[which(muts.ct$type=="ncpg")])
mtext(paste0("P=",signif(t$p.value,3)),3,cex=.7,line=-1)
dev.off()