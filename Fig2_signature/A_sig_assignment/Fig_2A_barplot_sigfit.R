#  read signatures
msi.sig<-read.table("Sigfit_sig_assign.txt",header = T,row.names = 1)
rownames(msi.sig)<-substring(rownames(msi.sig),1,12)

#  sample order
muts<-scan("MutS_sap",what="s")
mutl<-scan("MutL_sap",what="s")
sap.order<-c(intersect(muts,rownames(msi.sig)),intersect(mutl,rownames(msi.sig)))
n1<-length(intersect(muts,rownames(msi.sig)))
n2<-length(intersect(mutl,rownames(msi.sig)))
msi.sig<-msi.sig[sap.order,]
msi.sig<-subset(msi.sig,select = c(1,6,19,26,31,49)) #SBS1, SBS6, SBS15, SBS21, SBS26 and SBS44
msi.sig<-as.data.frame(msi.sig)
msi.sig$other<-(1-as.numeric(apply(msi.sig,1,sum)))

#  group
group<-data.frame(type=c(rep(1,n1),rep(2,n2)))
rownames(group)<-sap.order

#====================================
#      reorder
#===================================
# msi.sig<-msi.sig[rev(order(msi.sig[,1]+msi.sig[,3])),]
msi.sig<-msi.sig[rev(order(msi.sig[,1]+msi.sig[,2])),]
rname<-rownames(msi.sig)
group<-group[rname,]
cairo_pdf("Fig_2A_barplot_sigfit.pdf",9,3)
ly<-layout(matrix(c(1,2,3,3),ncol=2),heights = c(6,1),widths = c(20,2.5))

#============================
#      plot signatures
#============================
par(mar=c(0,1,.5,0),mgp=c(0.1,0.2,0))
cl<-brewer.pal(8, "Set2")[c(1,3:8)]
cl.six<-cl
tmp<-t(msi.sig)
barplot(tmp,col=cl,border=cl,space=0.2,axes=F,axisnames=F,ylab="Proportion",xpd=T)
at<-seq(0,1,0.2)
axis(side=2, at=at, labels=at, cex.axis=.8, xpd=TRUE,line=-0.9,tcl=-.2)

#=============================
#   plot mutants
#=============================
tmp<-as.matrix(group)
cl<-c("#ef8a62","#67a9cf")
cl.muts<-cl
par(mar=c(0.5,2.65,.5,1.75),mgp=c(2,0,0))
c<-image(tmp,col=cl,axes=F,font=2)
grid(dim(tmp)[1], dim(tmp)[2], col="white", lty=1, lwd=.5)
axis(at=seq(0,1,length.out=dim(tmp)[2]),labels="Mutant",side=2,las=3,tick=F,cex.axis=1,font=1,line=.4)

#=============================
#   plot legend
#=============================
##empty figure SIX TYPE
par(mgp=c(2,0,0),mar=c(0,0,0,.5))
plot(0,type="n",bty="n",ann=F,axes=F)
legend(x=0.45,y=.8,legend=c("SBS1","SBS6","SBS15","SBS21","SBS26","SBS44","Other"),x.intersp=0.3,y.intersp=1.2,fill=cl.six,horiz=F,xpd=T, box.col=NA,border=col,ncol=1)

# mutants
legend(x=0.45,y=-0.72,legend=c(paste0("MutS","\u03b1"),paste0("MutL","\u03b1")),y.intersp=1.2,cl.muts,horiz=F,xpd=T, box.col=NA,border=col,ncol=1,x.intersp=0.3,text.width =0.3)
dev.off()

