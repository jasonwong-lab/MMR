library(corrplot)
library(FactoMineR)
library(factoextra)
library(plotrix)
library(RColorBrewer)

#----------------------read data -------------------
all.data<-read.table("TCGA_MSI_anno.txt",header = T,row.names = 1)
all.data$muts<-as.character(all.data$muts)
all.data<-subset(all.data,muts!="Double")
order_name<-rownames(all.data)
group<-as.data.frame(all.data$muts)
rownames(group)<-order_name
colnames(group)<-"type"

#---------------------read trinucleotide------------
sig<-read.table("Tri_nucleotide.mat",row.names = 1,header = T,check.names = F)
rownames(sig)<-substring(rownames(sig),1,12)
sig<-as.data.frame(t(apply(sig,1,prop.table)))
sig.mat<-as.matrix(sig)
sig.mat<-sig.mat[order_name,]
sig.mat1<-sig.mat

#--------------------read 6 types of mutations
muttype<-read.table("Mut_type_six.txt",row.names = 1,header = T)
muttype<-subset(muttype,select=c(7:12))
rownames(muttype)<-substring(rownames(muttype),1,12)
muttype<-muttype[order_name,]

#---------PCA for 96 mut types of mutation ----------------
res.pca <- PCA(X = sig.mat1, scale.unit = T,graph = F)
dev.new()
pdf("Fig_1B.pdf")
par(mar=c(4,4,1,0.5),mgp=c(2,.5,0))
cl<-brewer.pal(6,"Set2")
tmp<-res.pca$ind$coord
data<-tmp[,1:2]
plot(data,col="white",xlab="PC1",ylab="PC2")
n<-nrow(data)
c.muts<-which(all.data$muts==1)
c.mutl<-which(all.data$muts==0)

for(i in c.mutl){
  floating.pie(data[i,1],data[i,2],as.numeric(muttype[i,]),radius=0.4,col=cl,lty = 1,lwd=.5)
}

for(i in c.muts){
  floating.pie(data[i,1],data[i,2],as.numeric(muttype[i,]),radius=0.4,col=cl,lty = 1,lwd=.1)
  draw.circle(data[i,1],data[i,2],radius=0.4,nv=100,border="#e41a1c",col=NA,lty=1,density=NULL, angle=45,lwd=2.5)
}

x<-max(data[,1])-2
y<-min(data[,2])

usr <- par("usr")
x<-usr[1]
y<-usr[3]*0.8

p1<-floating.pie(x*0.7,y,rep(1,6),radius=1,col=cl,lwd=.5)
pie.labels(x*0.7,y,p1,radius=1.1,labels=c("C>A","C>T","C>G","T>A","T>C","T>G"),cex=.8)
draw.circle(x*0.4,y,radius=0.7,nv=100,border=1,col=NA,lty=1,density=NULL, angle=45,lwd=1)
draw.circle(x*0.2,y,radius=0.7,nv=100,border="#e41a1c",col=NA,lty=1,density=NULL, angle=45,lwd=1)

dev.off()