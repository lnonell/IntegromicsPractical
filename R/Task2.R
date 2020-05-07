#Task 2: Alessia Squitieri 
#Purpose: Prepare expression data 
#input: Dataframe of function1 and cancer code 
#output: Dataframe with normalized (TMM) expression data 

##########################################################################
# Function task2:
##########################################################################

## Example of cancer type: 
cancer_type <- "TCGA-LUAD"

task2<-function(cancer_type, clinical.table){
  #libraries
  library(TCGAbiolinks)
  library(edgeR)
  library(SummarizedExperiment)
  library(biomaRt)
  
  #check arguments  
  stopifnot(is.character(cancer_type))
  stopifnot(is.data.frame(clinical.table))
  ########################################################################
  # 1. Download quey.rna:
  ########################################################################
  query.rna <- GDCquery(project = cancer_type, 
                      barcode = clinical.table$barcode,
                      data.category = "Transcriptome Profiling", 
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")
  getResults(query.rna) #data that will be downloaded
  GDCdownload(query.rna)
  TCGA.exp <- GDCprepare(query.rna) #this is creating the R object
  cat("Getting expression data...")
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

  #### Changing the colnames with the patients'barcode! 
  colnames(sub_data) <- substr(colnames(sub_data),1,12)
  
  #### Checking if there are some duplicates in the patients! If there are, the function will select the randomly 
  #### only the first patient! 
  if(any(duplicated(colnames(sub_data))) == TRUE){
  dupl <- which(duplicated(colnames(sub_data)))
  sub_data <- sub_data[,-dupl]
  }
  
  dim(sub_data)
  
  ########################################################################
  # 4. Normalization: 
  ########################################################################
  
  #### Now we have to normalize the data, we use the package edgeR and first of all
  #### it be created the DEGlist with the data that have to be normalized 
  cat("Your data are being normalized")
  sub_data <- DGEList(sub_data)
  
  #### Then is used the function calcNormFactor, by which the TMM normalizations is applied to our data! 
  norm_data <- calcNormFactors(sub_data)
  expression.table <- as.data.frame(norm_data$counts)
  cat("Done!")
  
  #### The final object is an expression table, with samples in columns and genes in rows! 
  return(expression.table)
  }

  ####################### To test the function ###########################
  LUAD.exp <- task2(cancer_type = "TCGA-LUAD" , clinical.table = LUAD.pts)
  KIRK.exp <- task2(cancer_type = "TCGA-KIRC", clinical.table = KIRK.pts)
  HNSC.exp <- task2(cancer_type = "TCGA-HNSC", clinical.table = HNSC.pts)
  STAD.exp <- task2(cancer_type = "TCGA-STAD", clinical.table = STAD.pts)
  LUSC.exp <- task2(cancer_type = "TCGA-LUSC", clinical.table = LUSC.pts)
  KICH.exp <- task2(cancer_type = "TCGA-KICH", clinical.table = KICH.pts)
  SKCM.exp  <- task2(cancer_type = "TCGA-SKCM", clinical.table = SKCM.pts)
  KIRP.exp <- task2(cancer_type = "TCGA-KIRP", clinical.table = KIRP.pts)
  ESCA.exp <- task2(cancer_type = "TCGA-ESCA", clinical.table = ESCA.pts)

