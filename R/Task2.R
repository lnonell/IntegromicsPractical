#Task 2: Alessia Squitieri 
#Purpose: prepare mRNA data
#input: dataframe of function1 
#output: dataframe with normalized expression data 

# PROOF function2
##########################################################################
library(TCGAbiolinks)
library(edgeR)
library(SummarizedExperiment)

task2<-function(){
  return()
}

## defining variable of interest from the dataframe of function1 : 
barcode <- LUAD.table$barcode
cancer_type <- "TCGA-LUAD"

# download quey.rna 
query.rna <- GDCquery(project = cancer_type, 
                      barcode = barcode,
                      data.category = "Transcriptome Profiling", 
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")

getResults(query.rna) #data that will be downloaded
GDCdownload(query.rna)

TCGA.exp <- GDCprepare(query.rna) #this is creating the R object

# subset sample
samples.exp <- substr(colnames(TCGA.exp),1,12)
length(samples.exp) # we obtain 23 samples instead of 20, some there are some replicates: 
length(unique(samples.exp))

s <- colnames(TCGA.exp)[!duplicated(substring(colnames(TCGA.exp), 1, 12))]

# matrix with cols = samples and rows = genes
data <- assay(TCGA.exp)

# to normalize data with edger we first created a DEGlist and then we used the function calcNormFactors:
data <- DGEList(data)
tmm <- calcNormFactors(data)
.
