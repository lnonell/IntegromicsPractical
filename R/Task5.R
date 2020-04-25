#Task 5 Óscar Baeza
#Purpose: Perform correlations (Spearman) (CN/mRNA and mRNA/methylation) and retrieve a table with good correlated genes and regression plots.
#Input: 2 data.frames with format being samples in columns and genes (HUGO Symbol) in rows (mRNA data),
# samples in columns and genes (HUGO Symbol) in rows containing the segment mean of each gene and sample (CN data),
# samples in columns and genes (HUGO Symbol) in rows containing the beta values mean of each gene and sample (methilation data).
#Output: data.frame with good correlations between data.frames (in absolute value >0.67)

#### FALTAN PLOTS Y PDF ####

task5<-function(df_mRNA, df_other, df_samples){
  
  #check arguments  
  stopifnot(is.data.frame(df_mRNA))
  stopifnot(is.data.frame(df_other))
  stopifnot(is.data.frame(df_samples))
  
  #...
  colnames(df_mRNA) <- gsub("-",".",df_samples$barcode)
  colnames(df_mRNA)==colnames(df_other)
  df_mRNA.s<-df_mRNA[,colnames(df_other)]
  (colnames(df_mRNA.s)==colnames(df_other))
  
  
  common <- intersect(sort(rownames(df_other)),sort(rownames(df_mRNA.s)))
  df_other.c <- df_other[common,] # give you common rows in df 1  
  df_mRNA.c <- df_mRNA.s[common,] # give you common rows in df 2
  
  
  resultsComb<-"/ResultsComb"
  dir.create(file.path(getwd(), resultsComb),showWarnings = FALSE)
  
  for (i in df_mRNA.c){
    #print(i)
    #miRNA.genes<-miRNAGenes(i)
    common.genes<-intersect(rownames(df_other.c),rownames(df_mRNA.c))
    #now correlations 
    lng<-length(common.genes)
    if (lng>0){
      cor.rho<-array(NA,lng)
      cor.pval<-array(NA,lng)
      #pdf(file.path(paste(getwd(),resultsComb,paste(mRNA,".corr.mRNA.miRNA.pdf",sep = ""),sep="/")))
      for (j in 1:lng){
        #miRNA.id<-rownames(df_other.c[df_other.c==j,])
        mRNA<-common.genes[j]
        mRNA.id<-rownames(df_mRNA.c[rownames(df_mRNA.c) == mRNA,])[1]
        y <- as.vector(unlist(df_other.c[mRNA.id,]))
        x <- as.vector(unlist(df_mRNA.c[mRNA.id,]))
        cor<-cor.test(x,y, method = "spearman",exact=FALSE)
        cor.pval[j]<-cor$p.value
        cor.rho[j]<-cor$estimate 
        #plot(x,y,main=mRNA,xlab="log2RMA expression",ylab="log2miRMA expression",type="p",xlim=c(0,16),ylim=c(0,16),cex=0.8)
        #fit <- lm(y ~ x)
        #abline(fit, col="chartreuse3",xlim=c(0,16))        
      }
      #dev.off()  #close pdf file
      cor.table<-data.frame(row.names = common.genes,"Rho"=as.vector(cor.rho),"pval"=as.vector(cor.pval))
      
      qqplot(x = cor.table$pval,y = cor.table$Rho, main = "P-value vs Rho plot")
      qqplot(x = cor.table$pval.adj,y = cor.table$Rho, main = "Adjusted P-value vs Rho plot")
      
      rownames(cor.table) <- common.genes
      pval.adj <- p.adjust(cor.table$pval,method = "fdr")
      cor.table$pval.adj <- as.vector(pval.adj)
      cor.sig.table <- subset(x = cor.table, abs(cor.table$Rho) >= 0.67 & cor.table$pval.adj < 0.05)
      
      
      #createSheet(wb, name = miRNA.id )
      #write.csv2(cor.table,file=file.path(resultsComb,paste("final_corr_table.csv",sep=".")))
      #writeWorksheet(wb,cor.table,sheet = miRNA.id, startRow = 1, startCol = 1, header=TRUE,rownames=NULL)
    }
  }  
  
  return(cor.sig.table)
}


