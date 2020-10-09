library(forestplot)
r.name<-c("H3K36me3","H3K27me3","H3K9me3","H3K4me3","H3K4me1","H3K27ac","DNase HS","Rep.timing","MBD4 ChIP")
c.name<-c("or","lower","upper","pval")

# mbd4
mydata<-read.table("data/logit_regr_or.txt",sep="\t")
mydata<-mydata[-1,]
rownames(mydata)<-r.name
colnames(mydata)<-c.name
clrs <- fpColors(box="royalblue",line="darkblue", summary="royalblue")
tabletext <-data.frame(c("", rownames(mydata)),c("HR", sprintf("%.2f", mydata[,"or"])),c("P-value",as.character(mydata$pval)))
tabletext<-as.matrix(tabletext)
ci<-mydata[,1:3]
#dev.new()
pdf("mbd4_forestplot.pdf",5,3)
forestplot(tabletext,rbind(rep(NA, 3),ci),col=clrs,xticks = c(0.8,0.9,1,1.1,1.2),xlog=TRUE,new_page = TRUE,hrzl_lines = list("2" = gpar(lty=1)),vertices = TRUE,xlab="Hazard Ratio",graph.pos = 3,txt_gp = fpTxtGp(label = list(gpar(fontfamily = ""),gpar(fontfamily = "",col = "#660000")),ticks = gpar(cex=1),xlab  = gpar(cex = 1)),title="MBD4 mutants")
dev.off()