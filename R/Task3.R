#Task 3: Ariadna Cilleros Portet  
#Purpose: Download, preapre and treat CNV data.
#input: Dataframe with the metadata realted to our samples (output task1) and TCGA cancer type (ex: TCGA-KIRC)
#output: Dataframe with the CNV data related to each gene and grouped by patient. 
#         The data has been converted to be more suitable: -1, 0, +1.

task3<- function(cancer, df_samples){

  #check arguments  
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
  df_ann <- df_ann[,c(2,3,4,6,11,12)] #take chr, start, end, segment_mean, HGNC and Patient columns
  
  ###########################
  ##
  ## STEP4: Group probes by gene and patient
  ##
  ###########################
  
  
  #group the genes by patients and if there is one gene not recorded for one patient, 
  #we will assume that it has two copies in the same, which sould be the normal situation: 
  
  for ( i in unique(df_ann$GeneSymbol)){
    gene.CN<-df_ann[df_ann$GeneSymbol==i,]
    gene.CN.byPat <- gene.CN %>% group_by(Patient) %>% summarise(CN = mean(Segment_Mean, na.rm = TRUE))
    colnames(gene.CN.byPat) <- c("Patient",i)
    
    #1. In case that at least one patient is missing
    if (dim(gene.CN.byPat)[1] != length(unique(df_ann$Patient))){ 
      patients <- setdiff(unique(df_ann$Patient), gene.CN.byPat$Patient) #variables with patients missing
      for (patient in patients){
        patient.miss <- c(patient, 0)
        gene.CN.byPat <- rbind(gene.CN.byPat, patient.miss) #add value for patient missing
      }
    }
    
    #2. Make sure that we have the same order for patients in the data frame
    gene.CN.byPat <- as.data.frame(gene.CN.byPat[match(gene.CN.byPat$Patient, unique(df_ann$Patient)),]) 
    
    #3. With the first gene we will create the dataframe to store all the results
    if (i == unique(df_ann$GeneSymbol)[1]){ 
      df_CN <- data.frame(gene.CN.byPat[2], row.names = gene.CN.byPat$Patient) 
    }
    
    #4. Add gene values for the rest of the genes
    else{
      df_CN <- cbind(df_CN, gene.CN.byPat[2]) 
    }
  }
  
  ###########################
  ##
  ## STEP5: Transform data
  ##
  ###########################
  
  #We will transform the values to -1, 0 and 1, being -1 less copies of the gene, 
  #0 the normal copies, and +1 more copies than the expected. 
  
  df_CN <- t(df_CN) #we will transpose to have genes in rows and patients in columns 
  
  for (i in 1:(dim(df_CN)[2])) df_CN[,i] <- as.integer(ifelse(2^(as.numeric(df_CN[,i])+1)>2.4,1,ifelse(2^(as.numeric(df_CN[,i])+1)<1.6,-1,0)))
  
  return(data.matrix(data.frame(df_CN, stringsAsFactors = FALSE)))
}
