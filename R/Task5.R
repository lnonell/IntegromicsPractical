#Task 5 Óscar Baeza
#Purpose: Perform correlations (Spearman) (CN/mRNA and mRNA/methylation) and retrieve a table with good correlated genes and regression plots.
#Input: 2 data.frames with format being samples in columns and genes (HUGO Symbol) in rows (mRNA data),
# samples in columns and genes (HUGO Symbol) in rows containing the segment mean of each gene and sample (CN data),
# samples in columns and genes (HUGO Symbol) in rows containing the beta values mean of each gene and sample (methilation data).
#Output: data.frame with good correlations between data.frames (in absolute value >0.67)

task5<-function(dataframe1, dataframe2, method){
  
  
  #check arguments  
  stopifnot(is.data.frame(dataframe1))
  stopifnot(is.data.frame(dataframe2))
  #...
  return(correlations)
}

df1 <- data.frame(matrix(rnorm(20), nrow=10, ncol = 10))
df2 <- data.frame(matrix(rnorm(30), nrow=10, ncol = 10))


cna_driven_genes <- find_cna_driven_gene( gene_cna=KIRP.cnv, gene_exp=KIRP.exp, gain_prop =0.67, loss_prop = 0.67,
                                          progress = TRUE, parallel = FALSE)

cor.test(df1, df2, method = "spearman")

correlations <- cor(t(df1), t(df2), method = "spearman")


cor.test(x = df1,y = df2,method = "spearman", ci = FALSE)

aa <- correlation(x = df1, y = df2, method = "spearman", cutoff = 0.67)


A <- as.matrix(df1)
B <- as.matrix(df2)
sapply(seq.int(dim(A)[1]), function(i) cor.test(A[1,], B[1,]))
unlist(a)

dm1 <- data.matrix(df1)
dm2 <- data.matrix(df2)
cc <- cor.test(dm1, dm2, method = "spearman")
cbind(cc$estimate,cc$p.value)



#### for sorting the results
targets <- correlations[correlations > 0.5 , ]

function(col2, col1) {
  cc <- cor.test(col2, col1, method = "spearman")
  cbind(cc$estimate,cc$p.value)
}


for (i in miRNAs){
  print(i)
  miRNA.genes<-miRNAGenes(i)
  miRNA.genes.deg<-intersect(miRNA.genes,mRNAs)
  #now correlations 
  lng<-length(miRNA.genes.deg)
  if (lng>0){
    cor.rho<-array(NA,lng)
    cor.pval<-array(NA,lng)
    miRNA.id<-rownames(res.miRNA[res.miRNA$miRNA==i,])
    y=as.vector(data.miRNA[miRNA.id,])
    for (j in 1:lng){
      mRNA<-miRNA.genes.deg[j]
      mRNA.id<-rownames(res.mRNA[!is.na(res.mRNA$Symbol) & res.mRNA$Symbol==mRNA,])[1]  
      x=as.vector(data.mRNA[mRNA.id,])
      cor<-cor.test(x,y, method = "spearman",exact=FALSE)
      cor.pval[j]<-cor$p.value
      cor.rho[j]<-cor$estimate 
      plot(x,y,main=mRNA,xlab="log2RMA expression",ylab="log2miRMA expression",type="p",xlim=c(0,16),ylim=c(0,16),col=cols,pch=pchs,cex=0.8)
      fit <- lm(y ~ x)
      abline(fit, col="chartreuse3",xlim=c(0,16))        
    }
    cor.table<-data.frame("miRNA ID"=rep(miRNA.id,lng),"miRNA"=rep(i,lng),miRNA.genes.deg,"Rho"=as.vector(cor.rho),"pval"=as.vector(cor.pval))
  }
}  
