mydata<-read.table("TCGA_MSI_anno.txt",header = T,row.names = 1,sep="\t")

#                   samples order                    
muts<-subset(mydata,select = 37)
muts$type<-rep(NA,nrow(muts))
muts$type[which(muts$muts=="")]=1
muts$type[which(muts$muts=="Double")]=2
muts$type[which(muts$muts==1)]=3
muts$type[which(muts$muts==0)]=4
order_name<-rownames(muts)[rev(order(muts$type))]

library(RColorBrewer)
dev.new()
cairo_pdf("Fig_1A.pdf",7,4.5)
ly<-layout(matrix(c(1:10,11,11,11,11),ncol=2),heights = c(6,6,4,4,1,1,1.3),widths = c(20,3))

#--------------------------------
#       SNV+INDEL
#--------------------------------
count<-subset(mydata,select = c(3,4))
count<-count[order_name,]
count$snv_log<-(count$total_mut)
count$indel_log<-(count$indel)
tmp<-t(as.matrix(count[,3:4]))
par(mar=c(0,1,1,0),mgp=c(0.1,0.2,0))
cl<-brewer.pal(3, "Paired")[1:2]
cl.snv<-cl
barplot(tmp,col=cl,border=cl,space=0,axes=F,axisnames=F,ylab="Mutation counts",xpd=T,log="y",xaxt="n")
a<-max(tmp[1,]+tmp[2,])
at<-seq(0,a,length.out=5)
axis(side=2, at=c(1,10,100,1000,10000), labels=c(1,"",100,"",10000), cex.axis=.8, xpd=TRUE,line=-1,tcl=-.2)

#--------------------------------
#     SIX TYPE MUTATION
#--------------------------------
six.mut<-subset(mydata,select = 5:10)
six.mut<-six.mut[order_name,]
tmp<-t(six.mut)
par(mar=c(0,1,.5,0),mgp=c(0.1,0.2,0))
cl<-brewer.pal(6, "Set2")
cl.six<-cl
barplot(tmp,col=cl,border=cl,space=0,axes=F,axisnames=F,ylab="Proportion",xpd=T)
at<-c(0,0.25,0.5,0.75,1)
axis(side=2, at=at, labels=at, cex.axis=.8, xpd=TRUE,line=-1,tcl=-.2)


#--------------------------------------
#                MMR MUTATION
#--------------------------------------
mmr.mut<-read.table("MMR_gene_mut_status.txt",header = T,row.names = 1)
mmr.mut<-mmr.mut[order_name,]
tmp<-as.matrix(mmr.mut)
genes<-substring(colnames(mmr.mut),1,4)
cl<-c("#d9d9d9",brewer.pal(4, "Set1"))
cl.mut<-cl
par(mar=c(0,2.65,.5,1.65),mgp=c(2,0,0))
c<-image(tmp,col=cl,axes=F,font=2)
grid(dim(tmp)[1], dim(tmp)[2], col="white", lty=1, lwd=.5)
axis(at=seq(0,1,length.out=dim(tmp)[2]),labels=genes,side=2,las=1,tick=F,cex.axis=.7,font=3,line=.4)


#--------------------------------------
#                MMR EXPRESSION
#--------------------------------------
mmr.exp<-subset(mydata,select = 31:34)
mmr.exp<-mmr.exp[order_name,]
tmp<-as.matrix(mmr.exp)
cl<-c("#377eb8","white","#e41a1c")
hcl<-colorRampPalette(cl)(100)
par(mar=c(0,2.65,.5,1.65),mgp=c(2,0,0))
min.max<-range(na.omit(tmp))
image(tmp,col=hcl,axes=F,font=2,zlim = c(-3,2.4))
axis(at=seq(0,1,length.out=dim(tmp)[2]),labels=genes,side=2,las=1,tick=F,cex.axis=.7,font=3,line=.4)


#--------------------------------------
#                MLH1 METHYL
#--------------------------------------
mlh1.met<-subset(mydata,select = 24)
mlh1.met<-as.data.frame(mlh1.met[order_name,])
rownames(mlh1.met)<-order_name
colnames(mlh1.met)<-"MLH1_me"
mlh1.met$type<-rep(NA,nrow(mlh1.met))
mlh1.met$type[which(is.na(mlh1.met$MLH1_me))]=0
mlh1.met$type[which(mlh1.met$MLH1_me>0.26)]=1
mlh1.met$type[which(mlh1.met$MLH1_me<0.26)]=2
tmp<-as.matrix(mlh1.met$type)
cl<-c("#d9d9d9","#d73027","#4575b4")
cl.met<-cl
par(mar=c(0,2.65,.5,1.65),mgp=c(2,0,0))
c<-image(tmp,col=cl,axes=F,font=2)
#grid(dim(tmp)[1], dim(tmp)[2], col="white", lty=1, lwd=.5)
axis(at=seq(0,1,length.out=dim(tmp)[2]),labels="MLH1",side=2,las=1,tick=F,cex.axis=.7,font=3,line=.4)


#--------------------------------------
#                TUMOR TYPE
#--------------------------------------
ttype<-subset(mydata,select = 1)
ttype<-as.data.frame(ttype[order_name,])
rownames(ttype)<-order_name
colnames(ttype)<-"cancer_type"
ttype$type<-rep(NA,nrow(ttype))
ttype$type[which(ttype$cancer_type=="UCEC")]=1
ttype$type[which(ttype$cancer_type=="COAD" |ttype$cancer_type=="READ" )]=2
ttype$type[which(ttype$cancer_type=="STAD")]=3
tmp<-as.matrix(ttype$type)
cl<-brewer.pal(5, "Set3")[3:5]
cl.ttype<-cl
par(mar=c(0,2.65,.5,1.65),mgp=c(2,0,0))
c<-image(tmp,col=cl,axes=F,font=2)
#grid(dim(tmp)[1], dim(tmp)[2], col="white", lty=1, lwd=.5)
axis(at=seq(0,1,length.out=dim(tmp)[2]),labels="Ttype",side=2,las=1,tick=F,cex.axis=.7,font=1,line=.4)


#--------------------------------------
#       MutL and MutS type
#--------------------------------------
muts<-subset(mydata,select = 37)
muts<-as.data.frame(muts[order_name,])
rownames(muts)<-order_name
muts$type<-rep(NA,nrow(muts))
muts$type[which(muts$muts=="")]=1
muts$type[which(muts$muts=="Double")]=2
muts$type[which(muts$muts==1)]=3
muts$type[which(muts$muts==0)]=4
tmp<-as.matrix(muts$type)
cl<-c("#d9d9d9","#7fc97f","#ef8a62","#67a9cf")
cl.muts<-cl
par(mar=c(0.5,2.65,.5,1.65),mgp=c(2,0,0))
c<-image(tmp,col=cl,axes=F,font=2)
#grid(dim(tmp)[1], dim(tmp)[2], col="white", lty=1, lwd=.5)
axis(at=seq(0,1,length.out=dim(tmp)[2]),labels="Mutant",side=2,las=1,tick=F,cex.axis=.7,font=1,line=.4)


#--------------------------------------
#           figure legend
#--------------------------------------
#empty figure SNV+INDEL
par(mgp=c(2,0,0),mar=c(0,0,0,1))
plot(0,type="n",bty="n",ann=F,axes=F)
legend(x=0.45,y=.5,legend=c("SNV","InDel"),x.intersp=0.3,y.intersp=2,fill=cl.snv,horiz=F,xpd=T, box.col=NA,border=col,ncol=1)

#empty figure SIX TYPE
par(mgp=c(2,0,0),mar=c(0,0,0,1))
plot(0,type="n",bty="n",ann=F,axes=F)
legend(x=0.45,y=.8,legend=c("C>A","C>T","C>G","T>A","T>C","T>G"),x.intersp=0.3,y.intersp=1,fill=cl.six,horiz=F,xpd=T, box.col=NA,border=col,ncol=1)

#empty figure MMR MUT
par(mgp=c(2,0,0),mar=c(0,0,0,1))
plot(0,type="n",bty="n",ann=F,axes=F)
legend(x=0.45,y=1,legend=c("Truncate","Missense","InFrame","Splice"),x.intersp=0.2,y.intersp=1,fill=cl.mut[2:5],horiz=F,xpd=T, box.col=NA,border=col,ncol=1,title="Mutation type:",title.adj=0.25,text.width = 1)

#empty others
par(mgp=c(2,0,0),mar=c(0,0,0,1))
plot(0,type="n",bty="n",ann=F,axes=F)
cl<-c("#d73027","#4575b4",brewer.pal(5, "Set3")[3:5],c("#d9d9d9","#7fc97f","#ef8a62","#67a9cf"))
lb<-c(paste0("\u03b2",">0.25"),paste0("\u03b2","<0.25"),"UCEC","CRC","STAD","ND","Both",paste0("MutS","\u03b1"),paste0("MutL","\u03b1"))

# MLH1 methylation
legend(x=0.45,y=1,legend=lb[1:2],y.intersp=0.8,fill=cl[1:2],horiz=F,xpd=T, box.col=NA,border=col,ncol=2,title="Methylation:",title.adj=0.3,x.intersp=0.1,text.width =0.3)

# Cancer type
legend(x=0.45,y=0.5,legend=lb[3:5],y.intersp=0.8,fill=cl[3:5],horiz=F,xpd=T, box.col=NA,border=col,ncol=2,title="Ttype:",title.adj=0.2,x.intersp=0.3,text.width =0.3)

# MutS and MutL mutants
legend(x=0.45,y=-0.1,legend=lb[6:9],y.intersp=0.8,fill=cl[6:9],horiz=F,xpd=T, box.col=NA,border=col,ncol=2,title="Mutants:",title.adj=0.2,x.intersp=0.3,text.width =0.3)
dev.off()