library(sigfit)
data("cosmic_signatures", package = "sigfit")
mydata<-read.table("data/tri_nucle_96.mat",header = T,row.names = 1,check.names = F)
colnames(mydata)<-colnames(cosmic_signatures)

sigfit::plot_spectrum(mydata,pdf_path = "cell_line_profile.pdf")

