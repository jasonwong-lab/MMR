dev.new()
cairo_pdf("Fig_2C-F_sig_boxplot.pdf",4,3.5)
layout(matrix(1:4,nrow=2,byrow = T))
par(mar=c(2,4,1,1),mgp=c(2,.5,0),bty="l")
cl<-c("#67a9cf","#ef8a62","#8da0cb")
ylim<-c(0,1)

#===============================================
#              TCGA MSI
#===============================================
#-------------------------read muts------------------------
msi.all<-read.table("data/TCGA_MSI_anno.txt",header = T,row.names = 1)
colnames(msi.all)[35]<-"sigB"
colnames(msi.all)[36]<-"sigA"
#-----------------------combine-----------------
all.data<-msi.all
all.data$muts<-as.character(all.data$muts)
all.data<-subset(all.data,muts!="Double")
#---------------------ALL--------------
mmr.all<-all.data
boxplot(sigA~muts,data=mmr.all,main="TCGA MMRd",ylab="SigA contribution",border="#cbd5e8",notch=F,names=c(paste0("MutL","\u03b1"),paste0("MutS","\u03b1")),ylim=ylim)
stripchart(sigA~muts,data=mmr.all, method = "jitter",pch = 20, add = TRUE, col =cl , cex = 1,vertical=T)
t<-t.test(mmr.all$sigA[which(mmr.all$muts==1)],mmr.all$sigA[which(mmr.all$muts==0)])
mtext(paste0("P=",signif(t$p.value,3)),3,cex=.7,line=-1)


#===============================================
#              MSK-CRC
#===============================================
#----------read muts mutl-------------
mydata<-read.table("data/MSK_CRC_anno.txt",row.names = 1,header = T)
msi.sap<-rownames(mydata)
mydata$muts<-rep(NA,nrow(mydata))
mydata$muts[which(mydata$mutl==1)]=0
mydata$muts[which(mydata$mutl==0)]=1
#-------------------read signature---------------
sig<-read.table("data/MSK_CRC_sig.txt",header = T,row.names = 1)
sig<-sig[msi.sap,]
colnames(sig)<-c("SigB_new","SigA_new")
mmr.all<-cbind(mydata,sig)
#--------------------sigA-------------------
boxplot(SigA_new~muts,data=mmr.all,main="MSK-CRC",ylab="SigA contribution",border="#cbd5e8",notch=F,names=c(paste0("MutL","\u03b1"),paste0("MutS","\u03b1")),ylim=ylim)
stripchart(SigA_new~muts,data=mmr.all, method = "jitter",pch = 20, add = TRUE, col =cl[1:2] , cex = 1,vertical=T)
t<-t.test(mmr.all$SigA_new[which(mmr.all$muts==1)],mmr.all$SigA_new[which(mmr.all$muts==0)])
mtext(paste0("P=",signif(t$p.value,3)),3,cex=.7,line=-1)


#===============================================
#              MSK-UCEC
#===============================================
#----------read muts mutl-------------
mydata<-read.table("data/MSK_UCEC_anno.txt",row.names = 1,header = T,sep="\t")
msi.sap<-rownames(mydata)
mydata$muts<-rep(NA,nrow(mydata))
mydata$muts[which(mydata$mutl==1)]=0
mydata$muts[which(mydata$mutl==0)]=1
#-------------------read signature---------------
sig<-read.table("data/MSK_UCEC_sig.txt",header = T,row.names = 1)
sig<-sig[msi.sap,]
colnames(sig)<-c("SigB_new","SigA_new")
mmr.all<-cbind(mydata,sig)
#--------------------sigA-------------------
boxplot(SigA_new~muts,data=mmr.all,main="MSK-UCEC",ylab="SigA contribution",border="#cbd5e8",notch=F,names=c(paste0("MutL","\u03b1"),paste0("MutS","\u03b1")),ylim=ylim)
stripchart(SigA_new~muts,data=mmr.all, method = "jitter",pch = 20, add = TRUE, col =cl[1:2] , cex = 1,vertical=T)
t<-t.test(mmr.all$SigA_new[which(mmr.all$muts==1)],mmr.all$SigA_new[which(mmr.all$muts==0)])
mtext(paste0("P=",signif(t$p.value,3)),3,cex=.7,line=-1)


#===============================================
#              DEPMAP
#===============================================
mmr.all<-read.table("data/DEPMAP_MSI_anno.txt",header = T,row.names = 1,sep = "\t")
mmr.all$muts<-as.character(mmr.all$muts)
mmr.all<-subset(mmr.all,muts!="Double")
#------------------------plot sigA---------------------
boxplot(sigA~muts,data=mmr.all,main="Depmap MMRd",ylab="SigA contribution",border="#cbd5e8",notch=F,names=c(paste0("MutL","\u03b1"),paste0("MutS","\u03b1")),ylim=ylim)
stripchart(sigA~muts,data=mmr.all, method = "jitter",pch = 20, add = TRUE, col =cl, cex = 1,vertical=T)
t<-t.test(mmr.all$sigA[which(mmr.all$muts==1)],mmr.all$sigA[which(mmr.all$muts==0)])
mtext(paste0("P=",signif(t$p.value,3)),3,cex=.7,line=-1)

dev.off()
