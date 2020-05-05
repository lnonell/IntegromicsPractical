#Task 3: Ariadna Cilleros Portet  
#Purpose: Download, preapre and treat CNV data.
#input: Dataframe with the metadata realted to our samples (output task1) and TCGA cancer type (ex: TCGA-KIRC)
#output: Dataframe with the CNV data related to each gene and grouped by patient. 
#         The data won't be converted by default (transform = FALSE), in case transform = TRUE, 
#         it will converted the data to be more suitable: -1, 0, +1.

task3<- function(cancer, df_samples, transform = FALSE){
  
  #check arguments  
  stopifnot(is.character(cancer))
  stopifnot(is.data.frame(df_samples))
  stopifnot(is.logical(transform))
  
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
  
  cat("\nDownloading CN data from ", cancer, " cancer dataset...\n")
  
  #prepare the data to be downloaded from TCGA
  query.CN <- GDCquery(project = cancer, 
                       data.category = "Copy Number Variation",
                       data.type = "Masked Copy Number Segment", 
                       barcode = df_samples$barcode)
  GDCdownload(query.CN)
  TCGA.CN<- GDCprepare(query.CN)
  
  #we will create a dataframe to help us to remove the duplicates
  df_barcode_sample <- data.frame(Barcode = unique(TCGA.CN$Sample), Patient = substr(unique(TCGA.CN$Sample),1,12) )
  
  if (any(duplicated(df_barcode_sample$Patient))== TRUE){
    dupl <- which(duplicated(df_barcode_sample$Patient))
    barcodes_dupl <- df_barcode_sample[c(dupl), ]$Barcode
  }
  
  #remove duplicated samples
  i = 1
  for (barcode in barcodes_dupl){
    if (i == 1){
      df_not_dupl <- TCGA.CN[TCGA.CN$Sample != barcode, ]
    }
    else{
      df_not_dupl <- df_not_dupl[df_not_dupl$Sample != barcode, ]
    }
    i = i+1
  }
  
  ###########################
  ##
  ## STEP2: Annotate hg38 genes
  ##
  ###########################
  
  cat("\nAnnotating h38 genes...")
  
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
  
  cat("\nFind the correspondence between the hg38 genes and our segments of CN...\n")
  
  #to make an easier overlap we need GRanges object from our annotation
  colnames(genes) <- c("Chr","Start","End", "GeneSymbol") 
  genes_GR <- makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
  
  #to make an easier overlap we need GRanges object from our CN probes 
  colnames(df_not_dupl) <- c("GDC_Aliquot","chr", "start", "end", "Num_Probes", "Segment_Mean", "Sample")
  df_GR <- makeGRangesFromDataFrame(df_not_dupl[,c(2,3,4)])#we only take chr, start and end to do the overlap
  
  #overlap the annotation with our CNV probes
  hits <- findOverlaps(genes_GR, df_GR, type="within") #hits found
  df_ann <- cbind(df_not_dupl[subjectHits(hits),],genes[queryHits(hits),]) 
  df_ann$Sample<-substr(df_ann$Sample,1,12) #create patient variable
  df_ann <- df_ann[,c(2,3,4,6,7,11)] #take chr, start, end, segment_mean, HGNC, patient and sample columns
  
  
  ###########################
  ##
  ## STEP4: Group genes by patient 
  ##
  ###########################
  
  #Problems: 
  #1) more than one measure per gene in one patient --> we will perfom the mean of the values
  #2) in the case that one patient has a missing value -->we will assume that it is 0 (no gains, no lost = two normal copies)
  
  pb <- txtProgressBar(min = 0, max = length(unique(df_ann$GeneSymbol)), style = 3)
  
  cat("\nGrouping genes by patient & treating missing values... Be patient, this might take a while...\n")
  
  for (gene in unique(df_ann$GeneSymbol)){
    
    #perform progress bar
    i<-i+1
    setTxtProgressBar(pb, i)
    
    gene.CN<-df_ann[df_ann$GeneSymbol==gene,]
    gene.CN.byPat <- data.frame()
    
    #1. Do the mean for more than one measure per gene in one sample
    gene.CN.byPat <- gene.CN %>% group_by(Sample) %>% summarise(CN = mean(Segment_Mean, na.rm = TRUE))
    colnames(gene.CN.byPat) <- c("Patient",gene)
    
    #2. Assune missing values in patients are normal situations 
    if (dim(gene.CN.byPat)[1] != length(unique(df_ann$Sample))){ 
      patients <- setdiff(unique(df_ann$Sample), gene.CN.byPat$Patient) #variables with patients missing
      for (patient in patients){
        patient.miss <- c(patient, 0)
        gene.CN.byPat <- rbind(gene.CN.byPat, patient.miss) #add value for patient missing
      }
    }
    
    #Make sure that we have the same order of patients
    gene.CN.byPat <- as.data.frame(gene.CN.byPat[match(gene.CN.byPat$Patient, unique(df_ann$Sample)),]) 
    
    #With the first gene we will create the dataframe to store all the results
    if (gene == unique(df_ann$GeneSymbol)[1]){ 
      df_CN <- data.frame(gene.CN.byPat[2], row.names = gene.CN.byPat$Patient) 
    }
    
    #Add gene values for the rest of the genes
    else{
      df_CN <- cbind(df_CN, gene.CN.byPat[2]) 
    }
  }
  
  ###########################
  ##
  ## STEP5: Preprocess final dataframe 
  ##
  ###########################
  
  df_CN <- t(df_CN) #wgenes rows, columns patients
  
  #Transofrm continues to categorical values: -1 (lost), 0 (2 normal copies), +1 (gain)
  
  if (transform == TRUE){
    cat("\nPerforming transformation...\n")
    for (i in 1:(dim(df_CN)[2])) df_CN[,i] <- as.integer(ifelse(2^(as.numeric(df_CN[,i])+1)>2.4,1,ifelse(2^(as.numeric(df_CN[,i])+1)<1.6,-1,0)))
  }
  
  df_final <- data.frame(df_CN, stringsAsFactors = FALSE)
  
  colnames(df_final) <- gsub(pattern = "\\.", replacement = "-", colnames(df_final))
  
  return(df_final) #return a data frame as output
  
}
