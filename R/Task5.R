#Task 5 Óscar Baeza
#Purpose: Perform correlations (Spearman) (CN/mRNA and mRNA/methylation) and retrieve a table with good correlated genes and regression plots.
#Input: 2 data.frames with format being samples in columns and genes (HUGO Symbol) in rows (mRNA data),
# samples in columns and genes (HUGO Symbol) in rows containing the segment mean of each gene and sample (CN data),
# samples in columns and genes (HUGO Symbol) in rows containing the beta values mean of each gene and sample (methilation data).
#Output: data.frame with good correlations between data.frames (in absolute value >0.67)


task5 <- function(df_mRNA, df_other, df_samples, num_plots, rhovalue_sig = 0.67, pvalue_sig = 0.05){
  
  # Check arguments  
  stopifnot(is.data.frame(df_mRNA))
  stopifnot(is.data.frame(df_other))
  stopifnot(is.data.frame(df_samples))
  stopifnot(is.numeric(num_plots))
  
  # Libary loading
  suppressPackageStartupMessages(library(TCGAbiolinks))
  suppressPackageStartupMessages(library(GenomicFeatures))
  suppressPackageStartupMessages(library(tidyverse))
  
  # We put all samples on same format
  colnames(df_mRNA) <- gsub("-",".",colnames(df_mRNA))
  
  # To check if all samples are in the same order and correct it
  if (all(colnames(df_mRNA)==colnames(df_other)) == FALSE) {
    df_mRNA.s<-df_mRNA[,colnames(df_other)]
    cat("\nThe samples have been reorganized...\n")
  }
  
  # We get all the common genes in both dataframes sorted
  common.genes <- intersect(sort(rownames(df_other)),sort(rownames(df_mRNA.s)))
  df_other.c <- df_other[common.genes,] # We select common rows in df other  
  df_mRNA.c <- df_mRNA.s[common.genes,] # We select common rows in df mRNA
  
  # We create the Results directory, hanging from the working directory
  cat("\nCreating /Results directory...\n")
  results_dir<-"/Results"
  dir.create(file.path(getwd(), results_dir),showWarnings = FALSE)
  
  cat("\nPerforming correlations, this can take a while...\n")
  
  
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
      
    }
    
  }
  
  # Assemble the correlation table
  cat("\nAssembling correlation table...\n")
  cor.table<-data.frame(row.names = common.genes,"Rho"=as.vector(cor.rho),"pval"=as.vector(cor.pval))
  rownames(cor.table) <- common.genes
  
  # We adjust the pvalues for multiple comparisons
  cat("\nAdjusting pvalues for multiple comparisons...\n")
  pval.adj <- p.adjust(cor.table$pval,method = "fdr")
  cor.table$pval.adj <- as.vector(pval.adj)
  
  # We make the Pvalue plots
  cat("\nCreating PDFs for the Pvalue plots...\n")
  pdf(file.path(paste(getwd(),results_dir,"pvalue.plots.pdf",sep="/")))
  qqplot(x = cor.table$pval,y = cor.table$Rho, main = "P-value vs Rho plot")
  qqplot(x = cor.table$pval.adj,y = cor.table$Rho, main = "Adjusted P-value vs Rho plot")
  dev.off()
  
  # We subset from the original correlation table to get significance
  cor.sig.table <- subset(x = cor.table, abs(cor.table$Rho) >= rhovalue_sig & cor.table$pval.adj < pvalue_sig)
  
  
  # We make the plots for only those that are significative and selecting the number by argument
  genes.sig <- rownames(cor.sig.table)
  
  if (length(genes.sig)>= num_plots) {
    lng <- length(genes.sig[1:num_plots])
    
    cat("\nCreating PNGs for the plots...\n")
    # Loop for doing png plots per gene
    for (j in 1:lng){
      gene<-genes.sig[j]
      gene.id<-rownames(KIRP.exp.c[rownames(KIRP.exp.c) == gene,])[1]
      y <- as.vector(unlist(KIRP.cnv.c[gene.id,]))
      x <- as.vector(unlist(KIRP.exp.c[gene.id,]))
      png(file.path(paste(getwd(),results_dir,paste(gene,".mRNA.CNV.correlations.pdf",sep = ""),sep="/")))
      plot(y,x,main=gene,xlab="Copy Number Variation",ylab="RNA expression",type="b")
      dev.off()
    }
  } else {
    cat("\nOnly",length(genes.sig),"plots can be generated, and you asked for",num_plots,"\n")
  }
  
  # Final function return
  cat("\n DONE \n")
  return(cor.sig.table)  
}  

