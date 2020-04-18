#Task 1 Monica Sanchez
#Purpose of the task: Prepare clinical data
#input: TCGA cancer code
#output: Dataframe with samples barcode, stageI or IV, years smoked

task1<-function(cancer_type){
  
    #load library
    library(TCGAbiolinks)

    ###
    #1 GET patients with rna, cn and met data
    ###
    
    if (class(cancer_type) != "character"){
      stop(print("The input 'cancer_type' must be the TCGA cancer code. E.g.: 'TCGA-LUAD'"))
    }
  
  
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
    
    common.patients <- intersect(substr(getResults(query.rna, cols = "cases"), 1, 12),
                                 intersect(substr(getResults(query.cn, cols = "cases"), 1, 12),
                                           substr(getResults(query.met, cols = "cases"), 1, 12)))
    
    
    ###
    #2 Get clinical data
    ###
    
    clinical.data <- GDCquery_clinic(project = cancer_type, type = "clinical")
    
    common.patients <- as.data.frame(common.patients)
    colnames(common.patients) <- "submitter_id"
    
    clinical.common.patients <- merge(clinical.data, common.patients, by = "submitter_id")
    
    ###
    #3 Select patients with tumoral stage I and IV
    ###
    
    #To get the column number of the columns of interest: patient_id, tumor stage and years smoked
    id_col <- grep("submitter_id", colnames(clinical.common.patients))
    ts_col <- grep("tumor_stage", colnames(clinical.common.patients))
    ys_col <- grep("^years_smoked", colnames(clinical.common.patients))
    
    #To get up to 10 patients with stage I and 10 patients with stage II
    stageI <- clinical.common.patients[grep("stage i$|stage ia|stage ib", clinical.common.patients$tumor_stage),
                                       c(id_col, ts_col, ys_col)]
    
    if (nrow(stageI) >= 10){
      stageI <- stageI[1:10,]
    }
    
    stageIV <- clinical.common.patients[grep("stage iv", clinical.common.patients$tumor_stage),c(id_col, ts_col, ys_col)]
    
    if (nrow(stageIV) >= 10){
      stageIV <- stageIV[1:10,]
    }
    
    #To change "stage ia" and "stage ib" to "stage i" (same for "stage iv")
    
    clinical.table <- rbind(stageI, stageIV)
    clinical.table$tumor_stage <- as.factor(clinical.table$tumor_stage)
    
    levels(clinical.table$tumor_stage)[levels(clinical.table$tumor_stage)=="stage ia"] <- "stage i"
    levels(clinical.table$tumor_stage)[levels(clinical.table$tumor_stage)=="stage ib"] <- "stage i"
    levels(clinical.table$tumor_stage)[levels(clinical.table$tumor_stage)=="stage iva"] <- "stage iv"
    levels(clinical.table$tumor_stage)[levels(clinical.table$tumor_stage)=="stage ivb"] <- "stage iv"
    
    #Prepare the final output
    
    rownames(clinical.table) <- c(1:nrow(clinical.table))
    colnames(clinical.table) <- c("barcode", "tumor_stage", "years_smoked")

return(clinical.table)
    
}



#########################
#Test the function
#########################

###
#0 INPUT: cancer_type
###

#To get all the project id's from TCGA
all_cancer_list <- TCGAbiolinks:::getGDCprojects()$project_id
cancer_list <- all_cancer_list[grep("ESCA|HNSC|KICH|KIRC|KIRP|LUAD|LUSC|SKCM|STAD", all_cancer_list)]
cancer_list

LUAD.table <- task1(cancer_type = "TCGA-LUAD")
KIRK.table <- task1(cancer_type = "TCGA-KIRC")
HNSC.table <- task1(cancer_type = "TCGA-HNSC")
STAD.table <- task1(cancer_type = "TCGA-STAD")
LUSC.table <- task1(cancer_type = "TCGA-LUSC")
KICH.table <- task1(cancer_type = "TCGA-KICH")
SKCM.table <- task1(cancer_type = "TCGA-SKCM")
KIRP.table <- task1(cancer_type = "TCGA-KIRP")
ESCA.table <- task1(cancer_type = "TCGA-ESCA")
