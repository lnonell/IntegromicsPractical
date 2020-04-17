#Task 2: Alessia Squitieri 
#Purpose: prepare mRNA data
#input: dataframe of function1 
#output: dataframe with normalized expression data 

##########################################################################
# PROOF function2
##########################################################################
library(TCGAbiolinks)
library(edgeR)
library(SummarizedExperiment)
library(biomaRt)
##########################################################################
cancer_type <- "TCGA-LUAD"
clinical.table$barcode <- barcode 

task2<-function(cancer_type, clinical.table){
  ########################################################################
  # download quey.rna 
  query.rna <- GDCquery(project = cancer_type, 
                      barcode = barcode,
                      data.category = "Transcriptome Profiling", 
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")
  getResults(query.rna) #data that will be downloaded
  GDCdownload(query.rna)
  TCGA.exp <- GDCprepare(query.rna) #this is creating the R object
  ########################################################################
  # Annotation 
  ensembl <- useMart("ensembl")
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
  annot_df <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), 
                  filters = "ensembl_gene_id", values = rownames(TCGA.exp), mart = mart)
  annot_df <- annot_df[order(annot_df[,2]),] #order by start position
  genes <- annot_df[order(annot_df[,1]),] #order by chromosome
  
  ########################################################################
  # Create matrix of interest and check for duplicates: 
  data <- assay(TCGA.exp)
  patient <- substr(colnames(data),1,12)
  
  if(any(duplicated(patient)) == TRUE){
  dupl <- which(duplicated(patient))
  data <- data[,-dupl]
  }
  ########################################################################
  # Normalization: 
  data <- DGEList(data)
  norm_data <- calcNormFactors(data)
  expression.table <- as.data.frame(norm_data$counts)
  return(expression.table)
}

