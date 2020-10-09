dev.new()
pdf("wgs_asym_ct_methyl.pdf",8,4)
layout(matrix(1:2,nrow=1,byrow = T))
par(mar=c(2,3,1,1),mgp=c(2,.5,0),bty="l")
ylim<-c(-1,1)

#============================================================================================
#               plot methylation and bias
#============================================================================================
p3.list<-scan("data/combine/msi_prime3_list",what="s")
p5.list<-scan("data/combine/msi_prime5_list",what="s")
bins<-scan("data/combine/bin_name",what="s")
n<-6
met3<-as.data.frame(matrix(rep(0,2*6),ncol=2))
for(i in 1:4){
  p3.tmp<-grep(paste0(bins[i],"_"),p3.list,value = T)
  p3.data<-read.table(p3.tmp,row.names = 1)
  p5.tmp<-grep(paste0(bins[i],"_"),p5.list,value = T)
  p5.data<-read.table(p5.tmp,row.names = 1)
  name<-intersect(rownames(p3.data),rownames(p5.data))
  rep.all<-cbind(p3.data[name,],p5.data[name,])
  rep.all$ctall<-rep.all[,11]+rep.all[,24]
  rep.all$gaall<-rep.all[,12]+rep.all[,23]
  final<-subset(rep.all,select=c("ctall","gaall"))
  met3<-final+met3
}

met4<-as.data.frame(matrix(rep(0,2*6),ncol=2))
for(i in 5:12){
  p3.tmp<-grep(paste0(bins[i],"_"),p3.list,value = T)
  p3.data<-read.table(p3.tmp,row.names = 1)
  p5.tmp<-grep(paste0(bins[i],"_"),p5.list,value = T)
  p5.data<-read.table(p5.tmp,row.names = 1)
  name<-intersect(rownames(p3.data),rownames(p5.data))
  rep.all<-cbind(p3.data[name,],p5.data[name,])
  rep.all$ctall<-rep.all[,11]+rep.all[,24]
  rep.all$gaall<-rep.all[,12]+rep.all[,23]
  final<-subset(rep.all,select=c("ctall","gaall"))
  met4<-final+met4
}

met3.data<-data.frame(ratio=log(met3$ctall/met3$gaall,2),type=rep("Low",6))
met4.data<-data.frame(ratio=log(met4$ctall/met4$gaall,2),type=rep("High",6))
all.data<-rbind(met3.data,met4.data)
boxplot(ratio~type,data=all.data,main="MMRd (WGS)",ylab="Log2(lagging/leading)",border="#cbd5e8",notch=F,ylim=ylim,names=c("Low methyl","High methyl"))
stripchart(ratio~type,data=all.data, method = "jitter",pch = 20, add = TRUE, col =c("#984ea3","#ff7f00") , cex = 1,vertical=T)
abline(h=0,lty=2,col="#bdbdbd")
t<-t.test(all.data$ratio[which(all.data$type=="Low")],all.data$ratio[which(all.data$type=="High")])
mtext(paste0("P=",signif(t$p.value,3)),3,cex=.7,line=-1)


#=====================POLE=============================
p3.list<-scan("data/combine/pole_prime3_list",what="s")
p5.list<-scan("data/combine/pole_prime5_list",what="s")
bins<-scan("data/combine/bin_name",what="s")
n<-9
sap<-c("TCGA-A6-6141","TCGA-AA-3555","TCGA-AA-3977","TCGA-AA-A00N","TCGA-AZ-4315","TCGA-CA-6717","TCGA-CA-6718","TCGA-EI-6917","TCGA-F5-6814")
met3<-as.data.frame(matrix(rep(0,2*n),ncol=2))
for(i in 1:4){
  p3.tmp<-grep(paste0(bins[i],"_"),p3.list,value = T)
  p3.data<-read.table(p3.tmp,row.names = 1)
  p5.tmp<-grep(paste0(bins[i],"_"),p5.list,value = T)
  p5.data<-read.table(p5.tmp,row.names = 1)
  name<-intersect(rownames(p3.data),rownames(p5.data))
  rep.all<-cbind(p3.data[name,],p5.data[name,])
  rep.all$ctall<-rep.all[,11]+rep.all[,24]
  rep.all$gaall<-rep.all[,12]+rep.all[,23]
  final<-subset(rep.all,select=c("ctall","gaall"))
  final<-final[sap,]
  met3<-final+met3
}
met4<-as.data.frame(matrix(rep(0,2*n),ncol=2))
for(i in 5:12){
  p3.tmp<-grep(paste0(bins[i],"_"),p3.list,value = T)
  p3.data<-read.table(p3.tmp,row.names = 1)
  p5.tmp<-grep(paste0(bins[i],"_"),p5.list,value = T)
  p5.data<-read.table(p5.tmp,row.names = 1)
  name<-intersect(rownames(p3.data),rownames(p5.data))
  rep.all<-cbind(p3.data[name,],p5.data[name,])
  rep.all$ctall<-rep.all[,11]+rep.all[,24]
  rep.all$gaall<-rep.all[,12]+rep.all[,23]
  final<-subset(rep.all,select=c("ctall","gaall"))
  final<-final[sap,]
  met4<-final+met4
}
met3.data<-data.frame(ratio=log(met3$ctall/met3$gaall,2),type=rep("Low",n))
met4.data<-data.frame(ratio=log(met4$ctall/met4$gaall,2),type=rep("High",n))

all.data<-rbind(met3.data,met4.data)
boxplot(ratio~type,data=all.data,main="POLE mutants (WGS)",ylab="Log2(lagging/leading)",border="#cbd5e8",notch=F,ylim=ylim,names=c("Low methyl","High methyl"))
stripchart(ratio~type,data=all.data, method = "jitter",pch = 20, add = TRUE, col =c("#984ea3","#ff7f00") , cex = 1,vertical=T)
abline(h=0,lty=2,col="#bdbdbd")
t<-t.test(all.data$ratio[which(all.data$type=="Low")],all.data$ratio[which(all.data$type=="High")])
mtext(paste0("P=",signif(t$p.value,3)),3,cex=.7,line=-1)

dev.off()
