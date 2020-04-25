#Task 8
#Purpose: Plot results of MFA, mixOmics, together with rawdata (mRNA, CN and methylation data) and significative correlation data.
#Filter previously each raw data set to get the 10% of the genes with the highest sd. The plot should be informative.
#input: data.frames obtained as output files in tasks 1 to 5. 100 most correlated variables obtained in tasks 6 and 7. Path to store Circos plot
#output: Circos plot in specified format. 


task8<-function(Cancer_table, mrna, cn, meth, cor.sig.table, mfa, mixomics, path){
  
  #check arguments  
  stopifnot(is.data.frame(mrna))
  stopifnot(is.data.frame(cn))
  stopifnot(is.data.frame(meth))
  stopifnot(is.data.frame(cor.sig.table))
  stopifnot(is.data.frame(mfa))
  stopifnot(is.data.frame(mixomics))
  
  # filter 10% genes by highest sd
  sd <- apply(mrna,1,sd)
  sd10 <- head(sort(sd,decreasing = TRUE), round(nrow(mrna)*0.1))
  mrna_final <- mrna[names(sd10),]
  
  sd <- apply(cn,1,sd)
  sd10 <- head(sort(sd,decreasing = TRUE), round(nrow(cn)*0.1))
  cn_final <- cn[names(sd10),]
  
  sd <- apply(meth,1,sd)
  sd10 <- head(sort(sd,decreasing = TRUE), round(nrow(meth)*0.1))
  meth_final <- meth[names(sd10),]
  
  
  
  return()
}