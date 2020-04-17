#Task 3 by Ariadna Cilleros 
#Purpose: Download, treat and prepare CNV data files from TCGA.
#input: Dataframe with the metadata of the samples/patients and TCGA project name (ex: TCGA-ACC)
#output: Dataframe with the CNV data summarized by gene and patient. The data has been converted
#         to have suitable values (-1, 0, 1) to interpret the CN alterations. 

task3<- function(cancer, df_samples){
  
  stopifnot(is.character(cancer))
  stopifnot(is.data.frame(df_samples))
  
  suppressPackageStartupMessages(library(TCGAbiolinks))
  suppressPackageStartupMessages(library(GenomicFeatures))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  suppressPackageStartupMessages(library(biomaRt))
  
  ###########################
  ##
  ## STEP1: Download TCGA data
  ##
  ###########################
  
  #prepare the data to be downloaded from TCGA
  barcode <- df_samples$barcode #barcode of the samples
  query.CN <- GDCquery(project = cancer, 
                       data.category = "Copy Number Variation",
                       data.type = "Masked Copy Number Segment", 
                       barcode = barcode)
  GDCdownload(query.CN)
  TCGA.CN<- GDCprepare(query.CN)
  
  ###########################
  ##
  ## STEP2: Annotate hg38 genes
  ##
  ###########################
  
  #Set up an gene annotation template to use
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genes <- genes(txdb)
  ensembl <- useMart("ensembl")
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  
  #obtain the information chr, start, end and HGNC name of all the genes annotated in hg38
  annot_df <- getBM(attributes = c("chromosome_name","start_position","end_position","hgnc_symbol"), 
                    filters = "entrezgene_id", values = genes$gene_id, mart = mart)
  
  annot_df<- annot_df[annot_df[,4]!="" & annot_df[,1] %in% c(1:22,"X","Y"),] #remove those not annotated 
  annot_df <- annot_df[order(annot_df[,2]),] #order by start position
  genes <- annot_df[order(annot_df[,1]),] #order by chromosome
  
  ###########################
  ##
  ## STEP3: Overlap probes with annotation by GR objects
  ##
  ###########################
  
  #to make an easier overlap we need GRanges object from our annotation
  colnames(genes) <- c("Chr","Start","End", "GeneSymbol") 
  genes_GR <- makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
  
  #to make an easier overlap we need GRanges object from our CN probes 
  colnames(TCGA.CN) <- c("GDC_Aliquot","chr", "start", "end", "Num_Probes", "Segment_Mean", "Sample")
  df_GR <- makeGRangesFromDataFrame(TCGA.CN[,c(2,3,4)])#we only take chr, start and end to do the overlap
  
  #overlap the annotation with our CNV probes
  hits <- findOverlaps(genes_GR, df_GR, type="within") #hits found
  df_ann <- cbind(TCGA.CN[subjectHits(hits),],genes[queryHits(hits),]) 
  df_ann$Patient<-substr(df_ann$Sample,1,12) #create patient variable
  df_ann <- df_ann[,c(2,3,4,6,11,12)] #take chr, start, end, segment_mean, HGNC and Patient
  
  ###########################
  ##
  ## STEP4: Group probes by gene and patient
  ##
  ###########################
  
  
  #we will run a for loop in which for every gene it will do a subset of all the probes corresponding to it and making 
  #a summary of it for each patient.
  df_CN <-  data.frame()
  j <- 0
  for ( i in unique(df_ann$GeneSymbol)){
    j <- j+1
    gene.CN<-df_ann[df_ann$GeneSymbol==i,]
    gene.CN.byPat <- gene.CN %>% group_by(Patient) %>% summarise(CN = mean(Segment_Mean, na.rm = TRUE))
    colnames(gene.CN.byPat) <- c("Patient",i)
    if (dim(gene.CN.byPat)[1] != length(unique(df_ann$Patient))){ 
      patients <- setdiff(gene.CN.byPat$Patient, unique(df_ann$Patient))
      for (patient in patients){
        patient.miss <- c(patient, 0.0000001)
        gene.CN.byPat <- rbind(gene.CN.byPat, patient.miss)
      }
    }
    gene.CN.byPat <- as.data.frame(gene.CN.byPat[match(gene.CN.byPat$Patient, unique(df_ann$Patient)),])
    if (i == unique(df_ann$GeneSymbol)[1]){ #in case of being the first gene running, we will create a data frame to have as rownames the Patients
      print(gene.CN.byPat)
      df_CN <- as.data.frame(gene.CN.byPat)
      df_CN <- column_to_rownames(.data = as.data.frame(df_CN), var = "Patient") %>% head()
    }
    else{
      print(gene.CN.byPat[2])
      df_CN <- cbind(df_CN, gene.CN.byPat[2]) #add gene values in the dataframe in case it is not the first gene running
    }
    cat(j, " gene of ", length(unique(df_ann$GeneSymbol)), "\n")
  }
  
  ###########################
  ##
  ## STEP5: Transform data
  ##
  ###########################
  
  #we will transform the values to -1, 0 and 1. 
  
  for (i in 1:(dim(df_CN)[2])) df_CN[,i] <- ifelse(2^(df_CN[,i]+1)>2.4,1,ifelse(2^(df_CN[,i]+1)<1.6,-1,0))
  
  return(df_CN)
}
