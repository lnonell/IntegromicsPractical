#Task 6: Blanca Rodriguez Fernandez 
#Purpose: Apply MFA to the mRNA, CN and methylation data comparing stage iv vs stage i.
#input: files of task1, task2, task3, task4 + path for output files
#outputs: dataframe 100 most correlated variables with PC1 and PC2 + plots 
task6 <- function(df_samples, df.rna, df.cn, df.met, pth = getwd()){
  
  if (pth == getwd()){
    warning(print("Default output file is your current working directory"))
  }
  ## check arguments
  stopifnot(is.data.frame(df_samples))
  stopifnot(is.data.frame(df.rna))
  stopifnot(is.data.frame(df.cn))
  stopifnot(is.data.frame(df.met))
  
  suppressPackageStartupMessages(library(FactoMineR))
  
  ## filter 10% genes by sd
  SD <- apply(df.rna,1,sd)
  top10sd <- head(sort(SD,decreasing = TRUE), round(nrow(df.rna)*0.1))
  rna.f <- df.rna[names(top10sd),]
  
  SD <- apply(df.cn,1,sd)
  top10sd <- head(sort(SD,decreasing = TRUE), round(nrow(df.cn)*0.1))
  cn.f <- df.cn[names(top10sd),]
  
  SD <- apply(df.met,1,sd)
  top10sd <- head(sort(SD,decreasing = TRUE), round(nrow(df.met)*0.1))
  met.f <- df.met[names(top10sd),]
  
  ## Data preparation (duda aqui)
  rna4MFA <- rna.f[!is.na(rna.f[,1]),]
  cn4MFA <- cn.f[!is.na(cn.f[,1]),]
  met4MFA <- met.f[!is.na(met.f[,1]),]
  
  ## stopifnot( colnames(rna4MFA), colnames(cn4MFA), colnames(met4MFA)) (duda aqui)
  ## Define conditions 
  cond <- as.factor(df_samples$tumor_stage)
  
  ## Length of the variables (genes)
  rna.l<-nrow(rna4MFA)
  cn.l<-nrow(cn4MFA)
  met.l<-nrow(met4MFA)
  
  # New data frame with individuals in rows and variables in columns
  data4Facto<-data.frame(as.factor(cond),t(rna4MFA),t(cn4MFA),t(met4MFA)) 
  rownames(data4Facto) <- paste(c(df_samples$barcode, cond, sep("_")))
  
  ## Apply MFA. # duda type of data 
  res.cond <- MFA(data4Facto, group=c(1,rna.l,cn.l,met.l), type=c("n","c","n","c"), 
                  ncp=5, name.group=c("cond","RNA","CN","MET"),num.group.sup=c(1), graph = FALSE) 
  pdf(paste(pth,"MFAplots.pdf",sep ="/"))
  plot(res.cond, choix = "ind")
  plot(res.cond, choix = "ind", partial="all")
  plot(res.cond, choix = "axes")
  dev.off()
  
  ## Return 100 most correlated variables with 
  ## PC1 (dimension 1 of global PCA)
  PC1 <- names(head(sort(res.cond$global.pca$var$cor[,1],decreasing = TRUE),100))
  ## PC2 (dimension 2 of global PCA)
  PC2 <- names(head(sort(res.cond$global.pca$var$cor[,2],decreasing = TRUE), 100))
  return(data.frame(PC1,PC2))
}
