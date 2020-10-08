library(sigfit)
#--------------get colnames----------------
data("variants_21breast", package = "sigfit")
counts_21breast <- sigfit::build_catalogues(variants_21breast)
mydata<-read.table("Tri_nucle.mat",header = T,row.names = 1)

#------------all msi--------------
sap<-scan("Sap_list",what="s")
mydata<-mydata[sap,]
colnames(mydata)<-colnames(counts_21breast)
mcmc_samples_extr <- sigfit::extract_signatures(counts = mydata,
                                                nsignatures = 2,
                                                iter = 5000, 
                                                seed = 1)
save(mcmc_samples_extr,file = "mcmc_samples_extr_sig2.RData")
extr_signatures <- sigfit::retrieve_pars(mcmc_samples_extr,
                                         par = "signatures")

sigfit::plot_spectrum(extr_signatures,pdf_path = "Denovo_sig.pdf")


