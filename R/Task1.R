#Task 1 Monica Sanchez
#Purpose of the task: Prepare clinical data
#input: TCGA cancer code (TCGA biolinks)
#output: Dataframe with samples, stageI or IV, years smoked

task1<-function(cancer_type){
  return(clinical.table)
}

library(TCGAbiolinks)

###
#1
###

#To get all the project id's from TCGA
all_cancer_list <- TCGAbiolinks:::getGDCprojects()$project_id

cancer_list <- cancer_list[grep("ESCA1|HNSC|KICH|KIRC|KIRP|LUAD|LUSC|SKCM|STAD", all_cancer_list)]

cancer_list[1]

#This is the input of the function
cancer_type <- "TCGA-LUAD"

#To check the datasets of each cancer type
ESCAdatasets <- TCGAbiolinks:::getProjectSummary(cancer_list[1]) #ESCA

ESCAdatasets$data_categories$data_category

###
#2
###

#Check that these samples have rna, cn and meth data:

query.rna <- GDCquery(project = cancer_type, 
                      data.category = "Transcriptome Profiling", 
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")
query.cn <- GDCquery(project = cancer_type, 
                     data.category = "Copy Number Variation", 
                     data.type = "Masked Copy Number Segment")
query.met <- GDCquery(project = cancer_type, 
                      data.category = "DNA Methylation", 
                      platform = "Illumina Human Methylation 450",
                      data.type = "Methylation Beta Value")


# Get all patients that have DNA methylation and gene expression.
common.patients <- intersect(substr(getResults(query.rna, cols = "cases"), 1, 12),
                            intersect(substr(getResults(query.cn, cols = "cases"), 1, 12),
                             substr(getResults(query.met, cols = "cases"), 1, 12)))

length(common.patients)#457

###
#3
###

#query the clinical data:
clinical.data <- GDCquery_clinic(project = cancer_type, type = "clinical")

#To get the column number of the columns of interest: patient_id, tumor stage and years smoked
id_col <- grep("submitter_id", colnames(clinical.data))
ts_col <- grep("tumor_stage", colnames(clinical.data))
ys_col <- grep("years_smoked", colnames(clinical.data))



clinical.common.patients <- merge(clinical.data, common.patients, by =clinical.data$submitter_id)


#To get 10 patients with stage I and 10 patients with stage II
stageI <- clinical.data[grep("stage ia|stage ib", clinical.data$tumor_stage),c(id_col, ts_col, ys_col)]
stageI <- stageI[1:10,] #implement an if to consider the case of less than 10 samples per group
stageIV <- clinical.data[grep("stage iv|stage iv", clinical.data$tumor_stage),c(id_col, ts_col, ys_col)]
stageIV <- stageIV[1:10,]

clinical.table <- rbind(stageI, stageIV)
clinical.table$tumor_stage <- as.factor(clinical.table$tumor_stage)

levels(clinical.table$tumor_stage)[levels(clinical.table$tumor_stage)=="stage ia"] <- "stage i"
levels(clinical.table$tumor_stage)[levels(clinical.table$tumor_stage)=="stage ib"] <- "stage i"


##################################################################

TCGA.chol.exp #class: RangedSummarizedExperiment 
colData(TCGA.chol.exp)
rowData(TCGA.chol.exp)
rowRanges(TCGA.chol.exp)
assays(TCGA.chol.exp)
head(assay(TCGA.chol.exp))


TCGA.chol.exp.mut <- GDCprepare(query,add.gistic2.mut = c("PTEN","FOXJ1"))
TCGA.chol.exp.mut #information has been added to previous results
colData(TCGA.chol.exp.mut)

#copy number 
# query <- GDCquery(project = "TCGA-CHOL", 
#                   data.category = "Copy Number Variation",
#                   data.type = "Copy Number Segment")
# GDCdownload(query)  

# query <- GDCquery(project = "TCGA-CHOL", 
#                   data.category = "Copy number variation",
#                   legacy=T,
#                   platform="Affymetrix SNP Array 6.0",file.type="hg19.seg")
# GDCdownload(query)  

query <- GDCquery(project = "TCGA-CHOL", 
                  data.category = "Copy Number Variation",
                  data.type = "Masked Copy Number Segment")
GDCdownload(query)  

TCGA.chol.CN <- GDCprepare(query) 

TCGA.chol.CN #tibble with alteraions by sample, chromosome, start and end
class(TCGA.chol.CN)
dim(TCGA.chol.CN)
head(TCGA.chol.CN)

#methylation data it tales long time
# query <- GDCquery(project = "TCGA-CHOL", 
#                   data.category = "DNA Methylation", 
#                   platform = "Illumina Human Methylation 450")
# GDCdownload(query)
# TCGA.chol.meth <- GDCprepare(query) 
# class(TCGA.chol.meth) #
# dim(TCGA.chol.meth)
# head(TCGA.chol.meth)

#clinical data, the following code does not work
# query <- GDCquery(project = "TCGA-CHOL",
#                   data.category = "Clinical")
# GDCdownload(query) 
# TCGA.chol.clin <- GDCprepare(query)  #it does not work

#instead use this
clinical <- GDCquery_clinic(project = "TCGA-CHOL", type = "clinical")

# Downloading and prepare using legacy
query <- GDCquery(project = "TCGA-GBM",
                  data.category = "Protein expression",
                  legacy = TRUE, 
                  barcode = c("TCGA-OX-A56R-01A-21-A44T-20","TCGA-08-0357-01A-21-1898-20"))
GDCdownload(query)
data <- GDCprepare(query, save = TRUE, 
                   save.filename = "gbmProteinExpression.rda",
                   remove.files.prepared = TRUE)
                   