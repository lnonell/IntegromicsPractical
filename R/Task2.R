#Task 2: Alessia Squitieri 
#Purpose: Prepare expression data 
#input: Dataframe of function1 and cancer code 
#output: Dataframe with normalized expression data 

##########################################################################
# Function task2:
##########################################################################

library(TCGAbiolinks)
library(edgeR)
library(SummarizedExperiment)
library(biomaRt)

## Example of cancer type: 
cancer_type <- "TCGA-LUAD"
## Define barcode variable before using the function, from the clinical.table obtained by function1:
barcode <- LUAD.table$barcode # clinical table obtained form function1 

task2<-function(cancer_type, clinical.table){
  ########################################################################
  # 1. Download quey.rna:
  ########################################################################
  query.rna <- GDCquery(project = cancer_type, 
                      barcode = barcode,
                      data.category = "Transcriptome Profiling", 
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")
  getResults(query.rna) #data that will be downloaded
  GDCdownload(query.rna)
  TCGA.exp <- GDCprepare(query.rna) #this is creating the R object
  
  ########################################################################
  # 2. Annotation HUGO:
  ########################################################################
  ensembl <- useMart("ensembl")
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
  annot <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), 
                  filters = "ensembl_gene_id", values = rownames(TCGA.exp), mart = mart)
  
  ## Create expression matrix
  data <- assay(TCGA.exp)
  
  ## Check order 
  data <- data[annot$ensembl_gene_id, ]
  
  ## Subset with the ones that have HGNC symbol
  sub_data <- data[annot[annot$hgnc_symbol != "",]$ensembl_gene_id, ]
  
  rownames(sub_data) <- annot[annot$hgnc_symbol != "",]$hgnc_symbol
  head(sub_data)

  ########################################################################
  # 3. Create matrix of interest and check for duplicates: 
  ########################################################################

  colnames(sub_data) <- substr(colnames(sub_data),1,12)
  
  if(any(duplicated(colnames(sub_data))) == TRUE){
  dupl <- which(duplicated(colnames(sub_data)))
  sub_data <- sub_data[,-dupl]
  }
  
  dim(sub_data)
  
  ########################################################################
  # 4. Normalization: 
  ########################################################################
  sub_data <- DGEList(sub_data)
  norm_data <- calcNormFactors(sub_data)
  expression.table <- as.data.frame(norm_data$counts)
  return(expression.table)
}

  ####################### To test the function ###########################
  LUAD.exp <- task2(cancer_type, clinical.table = LUAD.table)
  LUSC.exp <- task2(cancer_type = cancer_type, clinical.table = LUSC.table)
  SKCM.exp <- task2(cancer_type, clinical.table = SKCM.table)

