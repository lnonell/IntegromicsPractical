#Task 4: Cristina Tuni i Dominguez
#Purpose: Prepare methylation data to be downloaded and further treated.
#input: Data frame with the prepared clinical data from task 1, plus the cancer code from TCGA.
#output: Data frame with samples in columns and genes (HUGO Symbol) in rows containing the beta values mean of each gene and sample.

task4<-function(cancer, df_samples){
  
  stopifnot(is.character(cancer))
  stopifnot(is.data.frame(df_samples))
  
  #Libary loading
  
  suppressPackageStartupMessages(library(TCGAbiolinks))
  suppressPackageStartupMessages(library(GenomicFeatures))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  suppressPackageStartupMessages(library(biomaRt))
  
  #1. Downloading methylation data
  
  query.meth <- GDCquery(project = cancer, 
                         data.category = "DNA Methylation",
                         platform = "Illumina Human Methylation 450", 
                         barcode = df_samples$barcode)
  GDCdownload(query.meth)
  TCGA.meth<- GDCprepare(query.meth)
  
  return()
}