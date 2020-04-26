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
  
  metadataGenesInfo<-genes.info@elementMetadata
  
  #get sample information
  sample.info<- colData(TCGA.meth)
  
  #Probe, gene symbol, and position to TSS
  subGenesInfo<-metadataGenesInfo[,c(1:2,5)]
  
  identical(subGenesInfo[[1]],rownames(data)) #TRUE, they are ordered in the same way
  
  #binding of previous info and samples
  final<-cbind(subGenesInfo,data)
  
  
  #removeNA
  final<-final[complete.cases(final),]
  
  #remove those with "." in gene_symbol
  
  final<-final[!(final$Gene_Symbol=="."),]
  
  
  #phenodata<-TCGA.meth@elementMetadata
  
  #test<-TCGA.meth@rowRanges
  
  #testObj<-annotateGRanges(data,txdb) #test 
  
  ########################################################
  #2. Annotation
  ########################################################
  
  
  # Create template to use in getBM(values=...)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genes <- genes(txdb)
  ensembl <- useMart("ensembl")
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
  
  #obtain chr, start, end and HGNC name of all the genes annotated in hg38
  annot_df <- getBM(attributes = c("chromosome_name","start_position","end_position","hgnc_symbol"), 
                    filters = "entrezgene_id", values = genes$gene_id, mart = mart)
  
  
  return()
}
