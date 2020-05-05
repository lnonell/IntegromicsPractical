#Task 4: Cristina Tuni i Dominguez
#Purpose: Prepare methylation data to be downloaded and further treated.
#input: Data frame with the prepared clinical data from task 1, plus the cancer code from TCGA.
#output: Data frame with samples in columns and genes (HUGO Symbol) in rows containing the beta values mean of each gene and sample.

task4<-function(cancer, df_samples){
  
  stopifnot(is.character(cancer))
  stopifnot(is.data.frame(df_samples))
  
  ########################################################
  #Libary loading
  ########################################################
  
  suppressPackageStartupMessages(library(TCGAbiolinks))
  suppressPackageStartupMessages(library(SummarizedExperiment))
  suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  suppressPackageStartupMessages(library(GenomicFeatures))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(biomaRt))
  
  ########################################################
  #1. Downloading methylation data
  ########################################################
  
  query.meth <- GDCquery(project = cancer, 
                         data.category = "DNA Methylation",
                         platform = "Illumina Human Methylation 450", 
                         barcode = df_samples$barcode)
  GDCdownload(query.meth)
  TCGA.meth<- GDCprepare(query=query.meth, 
                         save=T, 
                         save.filename = "DNAmeth.rda", 
                         summarizedExperiment = T)
  
  # get expression matrix
  data<-assay(TCGA.meth)
  
  #get genes information
  genes.info<- rowRanges(TCGA.meth)
  
  #check if the data object and the genes.info object have the same order (they should)
  all(genes.info@ranges@NAMES==rownames(data)) #TRUE
  
  #from the previous GRanges, obtain useful dataframe
  rownames(genes.info)==rownames(data)
  
  chromosome<-genes.info@seqnames
  start_position<-genes.info@ranges@start
  
  gene_name<-genes.info@elementMetadata@listData$Gene_Symbol
  i=0
  for (i in 1:length(gene_name)){
    gene_name[i]<-strsplit(gene_name[i], "[;]")[[1]][1]
  }
  
  position_to_TSS<-genes.info@elementMetadata@listData$Position_to_TSS
  i=0
  for (i in 1:length(position_to_TSS)){
    position_to_TSS[i]<-strsplit(position_to_TSS[i], "[;]")[[1]][1]
  }
  
  
  genes.info_df<-data.frame(chromosome,start_position, gene_name, position_to_TSS, data)
  dim(genes.info_df) #485577 27
  head(genes.info_df)
  
  #remove those probes without annotated genes
  genes.info_df<-genes.info_df[genes.info_df$gene_name != ".", ]
  
  #transform the positions to numeric type
  genes.info_df$position_to_TSS<-as.numeric(as.character(genes.info_df$position_to_TSS))
  
  #aggregated data frame
  first_sample<-aggregate(TCGA.05.4390.01A.02D.1756.05 ~ gene_name, FUN = mean, data=genes.info_df)
  
  #remove the NA that are found across all columns
  genes.info_df<-genes.info_df[complete.cases(genes.info_df), ]
  
  #remove those probes without annotated genes
  genes.info_df<-genes.info_df[genes.info_df$gene_name != ".", ]
  
  #transform the positions to numeric type
  genes.info_df$position_to_TSS<-as.numeric(as.character(genes.info_df$position_to_TSS))
  
  #get sample information
  sample.info<- colData(TCGA.meth)
  
  ########################################################
  #2. Annotation
  ########################################################
 
  
  return()
}
