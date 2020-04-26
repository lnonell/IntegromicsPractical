#Task 5 Óscar Baeza
#Purpose: Perform correlations (Spearman) (CN/mRNA and mRNA/methylation) and retrieve a table with good correlated genes and regression plots.
#Input: 2 data.frames with format being samples in columns and genes (HUGO Symbol) in rows (mRNA data),
# samples in columns and genes (HUGO Symbol) in rows containing the segment mean of each gene and sample (CN data),
# samples in columns and genes (HUGO Symbol) in rows containing the beta values mean of each gene and sample (methilation data).
#Output: data.frame with good correlations between data.frames (in absolute value >0.67)

#### FALTAN PLOTS Y PDF ####

task5<-function(df_mRNA, df_other, df_samples, correlation, rhovalue_sig = 0.67, pvalue_sig = 0.05){
  
  # Check arguments  
  stopifnot(is.data.frame(df_mRNA))
  stopifnot(is.data.frame(df_other))
  stopifnot(is.data.frame(df_samples))
  
  # Libary loading
  suppressPackageStartupMessages(library(TCGAbiolinks))
  suppressPackageStartupMessages(library(GenomicFeatures))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  suppressPackageStartupMessages(library(biomaRt))
  
  # We put all samples on same format
  colnames(df_mRNA) <- gsub("-",".",colnames(df_mRNA))
  
  # To check if all samples are in the same order and correct it
  if (all(colnames(df_mRNA)==colnames(df_other)) == FALSE) {
    df_mRNA.s<-df_mRNA[,colnames(df_other)]
    print("The samples have been reorganized")
  }

  # We get all the common genes in both dataframes sorted
  common.genes <- intersect(sort(rownames(df_other)),sort(rownames(df_mRNA.s)))
  df_other.c <- df_other[comcommon.genes,] # We select common rows in df other  
  df_mRNA.c <- df_mRNA.s[common.genes,] # We select common rows in df mRNA
  
  # We create the Results directory, hanging from the working directory
  results_dir<-"/Results"
  dir.create(file.path(getwd(), results_dir),showWarnings = FALSE)
  
  for (i in df_mRNA.c){
    
    # Now we will perform the correlations 
    # Initialize the variables to store pvalues and rho coefficient
    lng <- length(common.genes)
    if (lng>0) {
      cor.rho<-array(NA,lng)
      cor.pval<-array(NA,lng)
      

      # Loop for doing correlations per gene, storing the rho and pvalue
      for (j in 1:lng){
        gene<-common.genes[j]
        gene.id<-rownames(df_mRNA.c[rownames(df_mRNA.c) == gene,])[1]
        y <- as.vector(unlist(df_other.c[gene.id,]))
        x <- as.vector(unlist(df_mRNA.c[gene.id,]))
        cor<-cor.test(x,y, method = "spearman",exact=FALSE)
        cor.pval[j]<-cor$p.value
        cor.rho[j]<-cor$estimate
        
        # # We open the PDF to store plots depending on the type of correlation
        # if (correlation == "mRNA.CN") {
        #   pdf(file.path(paste(getwd(),results_dir,paste(gene,".mRNA.CN.corr.pdf",sep = ""),sep="/")))
        #   plot(x,y,main=gene,xlab="log2RNA expression",ylab="Copy Number Variation",type="b",xlim=c(0,16),ylim=c(0,4),cex=0.8)
        #   # Close pdf file
        #   dev.off()
        # } else {
        #   pdf(file.path(paste(getwd(),results_dir,paste(gene,".mRNA.Methylation.corr.pdf",sep = ""),sep="/")))
        #   plot(x,y,main=gene,xlab="log2RNA expression",ylab="Methylation",type="b",xlim=c(0,16),ylim=c(0,16),cex=0.8)
        #   #fit <- lm(y ~ x)
        #   #abline(fit, col="chartreuse3",xlim=c(0,16))
        #   # Close pdf file
        #   dev.off()
          
        }
        
      }
      
      # Assemble the correlation table
      cor.table<-data.frame(row.names = common.genes,"Rho"=as.vector(cor.rho),"pval"=as.vector(cor.pval))
      rownames(cor.table) <- common.genes
    }
    
    # We adjust the pvalues for multiple comparisons
    pval.adj <- p.adjust(cor.table$pval,method = "fdr")
    cor.table$pval.adj <- as.vector(pval.adj)
    
    # 
    pdf(file.path(paste(getwd(),results_dir,"pvalue.plots.pdf",sep="/")))
    qqplot(x = cor.table$pval,y = cor.table$Rho, main = "P-value vs Rho plot")
    qqplot(x = cor.table$pval.adj,y = cor.table$Rho, main = "Adjusted P-value vs Rho plot")
    dev.off()
    
    cor.sig.table <- subset(x = cor.table, abs(cor.table$Rho) >= rhovalue_sig & cor.table$pval.adj < pvalue_sig)
    
  return(cor.sig.table)  
  }  
#}

