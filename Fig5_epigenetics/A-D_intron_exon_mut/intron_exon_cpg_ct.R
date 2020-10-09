exon.len<-10.439090
intron.len<-317.326271

pdf("intron_exon_cpg_ct.pdf",8,2)
layout(matrix(1:4,nrow = 1,byrow = T))
par(mar=c(4,4,4,1),mgp=c(2,.5,0))

#==========================cpg data==================
exp.exon<-scan("data/exp_exon.txt",what="s")
exp.exon<-as.numeric(exp.exon)/c(100,100,100,10,100,100)/exon.len

exp.intron<-scan("data/exp_intron.txt",what="s")
exp.intron<-as.numeric(exp.intron)/c(100,100,100,10,100,100)/intron.len

obs.exon<-scan("data/obs_exon.txt",what="s")
obs.exon<-as.numeric(obs.exon)/exon.len
obs.intron<-scan("data/obs_intron.txt",what="s")
obs.intron<-as.numeric(obs.intron)/intron.len

#----------------------------------------------
all.data<-matrix(c(obs.intron,exp.intron,obs.exon,exp.exon),ncol=4)
all.data<-signif(all.data,3)

sap<-scan("data/sap_name",what="s")
sap[3]<-"MBD4 mutants"
sap[4]<-"POLE mutants"
sap[5]<-"MMRd(MSI)"
sap[6]<-"MSS"
cl<-c("#fbb4ae","#b3cde3")

for (i in c(6,4,5,3)){
   a<-barplot(all.data[i,],main=paste0(sap[i]),ylab="CpG C>T/Mb",names.arg=c("Obs","Exp","Obs","Exp"),col=cl,border=cl,space=c(0.1,0.1,0.4,0.1))
  tmp<-all.data[i,]/8
text(x=a,y=all.data[i,]-tmp,labels = all.data[i,],cex=.8)
mtext("Intron",1,cex=.8,line=2,adj=.2)
mtext("Exon",1,cex=.8,line=2,adj=.8)
}
dev.off()
