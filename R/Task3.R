#Task 3: Ariadna Cilleros Portet  
#Purpose: Download, preapre and treat CNV data.
#input: Dataframe with the metadata realted to our samples (output task1) and TCGA cancer type (ex: TCGA-KIRC)
#output: Dataframe with the CNV data related to each gene and grouped by patient. 
#         The data will be converted by default (transform = TRUE) converted to be more suitable: -1, 0, +1.

task3<- function(cancer, df_samples, transform = TRUE){
  
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
  
  #prepare the data to be downloaded from TCGA
  query.CN <- GDCquery(project = cancer, 
                       data.category = "Copy Number Variation",
                       data.type = "Masked Copy Number Segment", 
                       barcode = df_samples$barcode)
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
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
  
  #obtain chr, start, end and HGNC name of all the genes annotated in hg38
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
  df_ann$Sample<-substr(df_ann$Sample,1,15) #change the barcode of the whole aliquote to the fraction of the sample
  df_ann <- df_ann[,c(2,3,4,6,7,11,12)] #take chr, start, end, segment_mean, HGNC, patient and sample columns
  
  
  #At this step we can face three problems: 1) having more than one segment that measures 
  #the same gene 2) having patients with more than one sample, and 3) having patients
  #with missing values for a fragment. For the first problem, we will allways get the 
  #first sample per patient, for the second we will perform the mean between the 
  #different values assigned to the same gene and to the same sample, and for the third 
  #problem, we will assign a 0, assuming that this patient hasn't have alteration on that gene. 
  
  for (gene in unique(df_ann$GeneSymbol)){
    gene.CN<-df_ann[df_ann$GeneSymbol==gene,]
    gene.CN.bySam <- data.frame()
    gene.CN.byPat <- data.frame()
    
    #1. Take the mean of the samples that can have more than one segment mean per gene 
    gene.CN.bySam <- gene.CN %>% group_by(Sample) %>% summarise(CN = mean(Segment_Mean, na.rm = TRUE))
    colnames(gene.CN.bySam) <- c("Sample",gene)
    gene.CN.bySam$Sample <- substr(gene.CN.bySam$Sample,1,12)
    
    #2. Take only one sample per patient in case of having duplicates
    for (patient in unique(gene.CN.bySam$Sample)){
      gene.CN.byPat <- rbind(gene.CN.byPat, gene.CN[gene.CN$Patient == patient, c(7,4)][1,])
    }
    colnames(gene.CN.byPat) <- c("Patient",gene)
    
    #3) In case that at least one patient is missing
    if (dim(gene.CN.byPat)[1] != length(unique(df_ann$Patient))){ 
      patients <- setdiff(unique(df_ann$Patient), gene.CN.byPat$Patient) #variables with patients missing
      for (patient in patients){
        patient.miss <- c(patient, 0)
        gene.CN.byPat <- rbind(gene.CN.byPat, patient.miss) #add value for patient missing
      }
    }
    
    #Make sure that we have the same order for patients in the data frame
    gene.CN.byPat <- as.data.frame(gene.CN.byPat[match(gene.CN.byPat$Patient, unique(df_ann$Patient)),]) 
    
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
  
  df_CN <- t(df_CN) #we will transpose to have genes in rows and patients in columns 
  
  #We will transform the values to -1, 0 and 1, being -1 less copies of the gene, 
  #0 the normal copies, and +1 more copies than the expected. 
  if (transform == TRUE){
    for (i in 1:(dim(df_CN)[2])) df_CN[,i] <- as.integer(ifelse(2^(as.numeric(df_CN[,i])+1)>2.4,1,ifelse(2^(as.numeric(df_CN[,i])+1)<1.6,-1,0)))
  }
  
  return(data.matrix(data.frame(df_CN, stringsAsFactors = FALSE))) #let's make sure that we get a numeric matrix
  
}
