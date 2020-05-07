#Task 5 Óscar Baeza
#Purpose: Perform correlations (Spearman) (CN/mRNA and mRNA/methylation) and retrieve a table with good correlated genes and regression plots.
#Input: 2 data.frames with format being samples in columns and genes (HUGO Symbol) in rows (mRNA data),
# samples in columns and genes (HUGO Symbol) in rows containing the segment mean of each gene and sample (CN data),
# samples in columns and genes (HUGO Symbol) in rows containing the beta values mean of each gene and sample (methilation data).
#Output: data.frame with good correlations between data.frames (in absolute value >0.67)



#########################
### CNV Correlations ###
#########################


task5.cnv <- function(df_mRNA, df_CNV, df_samples, num_plots = 4, rhovalue_sig = 0.67, pvalue_sig = 0.05){
  
  # Check arguments  
  stopifnot(is.data.frame(df_mRNA))
  stopifnot(is.data.frame(df_CNV))
  stopifnot(is.data.frame(df_samples))
  stopifnot(is.numeric(num_plots))
  
  # Libary loading
  suppressPackageStartupMessages(library(GenomicFeatures))
  suppressPackageStartupMessages(library(tidyverse))
  
  # We put all samples on same format
  colnames(df_mRNA) <- gsub("-",".",colnames(df_mRNA))
  colnames(df_CNV) <- gsub("-",".",colnames(df_CNV))
  
  # We remove any NaN data from both dataframes
  df_CNV <- df_CNV[complete.cases(df_CNV), ]
  df_mRNA <- df_mRNA[complete.cases(df_mRNA), ]
  
  # To check if all samples are in the same order and correct it
  if (all(colnames(df_mRNA)==colnames(df_CNV)) == FALSE) {
    df_mRNA<-df_mRNA[,colnames(df_CNV)]
    cat("\nThe samples have been reorganized...\n")
  }

  
  # We get all the common genes in both dataframes sorted
  common.genes <- intersect(sort(rownames(df_CNV)),sort(rownames(df_mRNA)))
  df_CNV.c <- df_CNV[common.genes,] # We select common rows in df other  
  df_mRNA.c <- df_mRNA[common.genes,] # We select common rows in df mRNA
  
  # We create the Results directory, hanging from the working directory
  cat("\nCreating Results_CNV directory...\n")
  results_dir<-"/Results_CNV"
  dir.create(file.path(getwd(), results_dir),showWarnings = FALSE)
  
  cat("\nPerforming correlations. This might take a while...\n")
  
  
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
      y <- as.numeric(unlist(df_CNV.c[gene.id,]))
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
  pdf(file.path(paste(getwd(),results_dir,"pvalue_plots.pdf",sep="/")))
  qqplot(x = cor.table$pval,y = cor.table$Rho, main = "P-value vs Rho plot")
  qqplot(x = cor.table$pval.adj,y = cor.table$Rho, main = "Adjusted P-value vs Rho plot")
  dev.off()
  
  # We subset from the original correlation table to get significance
  cor.sig.table <- subset(x = cor.table, abs(cor.table$Rho) >= rhovalue_sig & cor.table$pval.adj < pvalue_sig)
  
  # We order this table so we have the better adj.p value on top
  cor.sig.table <- cor.sig.table[order(cor.sig.table$pval.adj),]
  
  # We make the plots for only those that are significative and selecting the number by argument
  genes.sig <- rownames(cor.sig.table)
  
  if (length(genes.sig)>= num_plots) {
    lng <- length(genes.sig[1:num_plots])
    
    cat("\nCreating PNGs for the plots...\n")
    # Loop for doing png plots per gene
    for (j in 1:lng){
      gene<-genes.sig[j]
      gene.id<-rownames(df_mRNA.c[rownames(df_mRNA.c) == gene,])[1]
      y <- as.vector(unlist(df_CNV.c[gene.id,]))
      x <- as.vector(unlist(df_mRNA.c[gene.id,]))
      png(file.path(paste(getwd(),results_dir,paste(gene,"_mRNA_CNV_correlations.png",sep = ""),sep="/")))
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

####################### To test the function ###########################

LUAD.cnv.corr <- task5.cnv(df_mRNA = LUAD.exp, df_CNV = LUAD.cnv, df_samples = LUAD.pts)
KIRC.cnv.corr <- task5.cnv(df_mRNA = KIRK.exp, df_CNV = KIRC.cnv, df_samples = KIRC.pts)
HNSC.cnv.corr <- task5.cnv(df_mRNA = HNSC.exp, df_CNV = HNSC.cnv, df_samples = HNSC.pts)
STAD.cnv.corr <- task5.cnv(df_mRNA = STAD.exp, df_CNV = STAD.cnv, df_samples = STAD.pts)
LUSC.cnv.corr <- task5.cnv(df_mRNA = LUSC.exp, df_CNV = LUSC.cnv, df_samples = LUSC.pts)
KICH.cnv.corr <- task5.cnv(df_mRNA = KICH.exp, df_CNV = KICH.cnv, df_samples = KICH.pts)
SKCM.cnv.corr <- task5.cnv(df_mRNA = SKCM.exp, df_CNV = SKCM.cnv, df_samples = SKCM.pts)
KIRP.cnv.corr <- task5.cnv(df_mRNA = KIRP.exp, df_CNV = KIRP.cnv, df_samples = KIRP.pts)
ESCA.cnv.corr <- task5.cnv(df_mRNA = ESCA.exp, df_CNV = ESCA.cnv, df_samples = ESCA.pts)


################################
### Methilation Correlations ###
################################



task5.meth <- function(df_mRNA, df_meth, df_samples, num_plots = 4, rhovalue_sig = 0.67, pvalue_sig = 0.05){
  
  # Check arguments  
  stopifnot(is.data.frame(df_mRNA))
  stopifnot(is.data.frame(df_meth))
  stopifnot(is.data.frame(df_samples))
  stopifnot(is.numeric(num_plots))
  
  # Libary loading
  suppressPackageStartupMessages(library(GenomicFeatures))
  suppressPackageStartupMessages(library(tidyverse))
  
  # We put all samples on same format
  colnames(df_mRNA) <- gsub("-",".",colnames(df_mRNA))
  colnames(df_meth) <- gsub("-",".",colnames(df_meth))
  
  # We erase the totalmean column from the data.frame if it exists
  if ("totalmean" %in% colnames(df_meth) == TRUE) {
    df_meth <- subset(df_meth, select=-c(totalmean))
    cat("\nErasing totalmean column...\n")
  }
  
  # We remove any NaN data from both dataframes
  df_meth <- df_meth[complete.cases(df_meth), ]
  df_mRNA <- df_mRNA[complete.cases(df_mRNA), ]
  
  # To check if all samples are in the same order and correct it
  if (all(colnames(df_mRNA)==colnames(df_meth)) == FALSE) {
    df_mRNA<-df_mRNA[colnames(df_meth),]
    cat("\nThe samples have been reorganized...\n")
  }
  
  
  # We get all the common genes in both dataframes sorted
  common.genes <- intersect(sort(rownames(df_meth)),sort(rownames(df_mRNA)))
  df_meth.c <- df_meth[common.genes,] # We select common rows in df other  
  df_mRNA.c <- df_mRNA[common.genes,] # We select common rows in df mRNA
  

  # We create the Results directory, hanging from the working directory
  cat("\nCreating /Results_Meth directory...\n")
  results_dir<-"/Results_Meth"
  dir.create(file.path(getwd(), results_dir),showWarnings = FALSE)
  
  cat("\nPerforming correlations. This might take a while...\n")
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
      y <- as.vector(unlist(df_meth.c[gene.id,]))
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
  
  # We set colors to best see signifcant values
  cor.table$color <- "black"
  cor.table$color[cor.table$pval.adj < pvalue_sig  & abs(cor.table$Rho) >= rhovalue_sig]="red"
  #cor.table$color[abs(cor.table$Rho) >= rhovalue_sig]="blue"
  
  # We make the Pvalue plots
  cat("\nCreating PDFs for the Pvalue plots...\n")
  pdf(file.path(paste(getwd(),results_dir,"pvalue_plots.pdf",sep="/")))
  qqplot(x = cor.table$pval,y = cor.table$Rho, main = "P-value vs Rho plot",col = cor.table$color)
  qqplot(x = cor.table$pval.adj,y = cor.table$Rho, main = "Adjusted P-value vs Rho plot",col = cor.table$color)
  dev.off()
  
  # We subset from the original correlation table to get significance
  cor.table <- cor.table[,1:3]
  cor.sig.table <- subset(x = cor.table, abs(cor.table$Rho) >= rhovalue_sig & cor.table$pval.adj < pvalue_sig)
  
  # We order this table so we have the better adj.p value on top
  cor.sig.table <- cor.sig.table[order(cor.sig.table$pval.adj),]
  
  # We make the plots for only those that are significative and selecting the number by argument
  genes.sig <- rownames(cor.sig.table)
  
  if (length(genes.sig)>= num_plots) {
    lng <- length(genes.sig[1:num_plots])
    
    cat("\nCreating PNGs for the plots...\n")
    # Loop for doing png plots per gene
    for (j in 1:lng){
      gene<-genes.sig[j]
      gene.id<-rownames(df_mRNA.c[rownames(df_mRNA.c) == gene,])[1]
      y <- as.vector(unlist(df_meth.c[gene.id,]))
      x <- as.vector(unlist(df_mRNA.c[gene.id,]))
      png(file.path(paste(getwd(),results_dir,paste(gene,"_mRNA_Meth_correlations.png",sep = ""),sep="/")))
      plot(y,x,main=gene,xlab="Methilation Beta values",ylab="RNA expression",type="b")
      dev.off()
    }
  } else {
    cat("\nOnly",length(genes.sig),"plots can be generated, and you asked for",num_plots,"\n")
  }
  # Final function return
  cat("\n DONE \n") 
  return(cor.sig.table)
}


####################### To test the function ###########################

LUAD.meth.corr <- task5.meth(df_mRNA = LUAD.exp, df_meth = LUAD.meth, df_samples = LUAD.pts)
KIRC.meth.corr <- task5.meth(df_mRNA = KIRK.exp, df_meth = KIRC.meth, df_samples = KIRC.pts)
HNSC.meth.corr <- task5.meth(df_mRNA = HNSC.exp, df_meth = HNSC.meth, df_samples = HNSC.pts)
STAD.meth.corr <- task5.meth(df_mRNA = STAD.exp, df_meth = STAD.meth, df_samples = STAD.pts)
LUSC.meth.corr <- task5.meth(df_mRNA = LUSC.exp, df_meth = LUSC.meth, df_samples = LUSC.pts)
KICH.meth.corr <- task5.meth(df_mRNA = KICH.exp, df_meth = KICH.meth, df_samples = KICH.pts)
SKCM.meth.corr <- task5.meth(df_mRNA = SKCM.exp, df_meth = SKCM.meth, df_samples = SKCM.pts)
KIRP.meth.corr <- task5.meth(df_mRNA = KIRP.exp, df_meth = KIRP.meth, df_samples = KIRP.pts)
ESCA.meth.corr <- task5.meth(df_mRNA = ESCA.exp, df_meth = ESCA.meth, df_samples = ESCA.pts)

