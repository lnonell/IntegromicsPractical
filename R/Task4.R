#Task 4: Cristina Tuni i Dominguez
#Purpose: Prepare methylation data to be downloaded and further treated.
#input: Data frame with the prepared clinical data from task 1, plus the cancer code from TCGA.
#output: Data frame with samples in columns and genes (HUGO Symbol) in rows containing the beta values mean of each gene and sample.

task4<-function(cancer, df_samples){
  
  stopifnot(is.character(cancer))
  stopifnot(is.data.frame(df_samples))
  
  ########################################################
  #Libary loading
  ########################################################
  
  suppressPackageStartupMessages(library(TCGAbiolinks))
  suppressPackageStartupMessages(library(SummarizedExperiment))
  suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  suppressPackageStartupMessages(library(GenomicFeatures))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(biomaRt))
  suppressPackageStartupMessages(library(dplyr))
  
  
  ########################################################
  #1. Downloading methylation data
  ########################################################
  
  cat("\nDownloading methylation data from ",cancer, "data. This may take a while...\n")
  
  query.meth <- GDCquery(project = cancer, 
                         data.category = "DNA Methylation",
                         platform = "Illumina Human Methylation 450", 
                         barcode = df_samples$barcode)
  GDCdownload(query.meth)
  TCGA.meth<- GDCprepare(query=query.meth, 
                         save=T, 
                         save.filename = "DNAmeth.rda", 
                         summarizedExperiment = T)
  
  #we will create a dataframe to help us to remove the duplicates
  df_barcode_sample <- data.frame(Barcode = unique(TCGA.meth$sample), Patient = substr(unique(TCGA.meth$sample),1,12) )
  
  if (any(duplicated(df_barcode_sample$Patient))== TRUE){
    dupl <- which(duplicated(df_barcode_sample$Patient))
    barcodes_dupl <- df_barcode_sample[c(dupl), ]$Barcode
  }
  
  #remove duplicated samples
  i = 1
  for (barcode in barcodes_dupl){
    if (i == 1){
      df_not_dupl <- TCGA.meth[TCGA.meth$sample != barcode, ]
    }
    else{
      df_not_dupl <- df_not_dupl[df_not_dupl$sample != barcode, ]
    }
    i = i+1
  }
  
  # get expression matrix
  data<-assay(df_not_dupl)
  
  #get genes information
  genes.info<- rowRanges(df_not_dupl)
  
  ########################################################
  #2. Annotation
  ########################################################
  
  cat("\nDownloading annotation of human genes.\n")
  
  #Set up an gene annotation template to use
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genes <- genes(txdb)
  ensembl <- useMart("ensembl")
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
  
  #obtain chr, start, and HGNC name of all the genes annotated in hg38
  annot_df <- getBM(attributes = c("chromosome_name","start_position","end_position","hgnc_symbol"), 
                    filters = "entrezgene_id", values = genes$gene_id, mart = mart)
  
  annot_df<- annot_df[annot_df[,4]!="" & annot_df[,1] %in% c(1:22,"X","Y"),] #remove those not annotated 
  annot_df <- annot_df[order(annot_df[,2]),] #order by start position
  genes <- annot_df[order(annot_df[,1]),] #order by chromosome
  
  ########################################################
  #3. Data processing
  ########################################################
  
  cat("\nCreating a new object to work with.\n")
  
  #check if the data object and the genes.info object have the same order (they should)
  #all(genes.info@ranges@NAMES==rownames(data)) #TRUE
  
  #from the previous GRanges, obtain useful dataframe
  
  chr<-genes.info@seqnames
  
  start<-as.double(genes.info@ranges@start)
  end<-as.double(start+1)
  
  gene_name<-genes.info@elementMetadata@listData$Gene_Symbol
  i=0
  for (i in 1:length(gene_name)){
    gene_name[i]<-strsplit(gene_name[i], "[;]")[[1]][1]
  }
  
  
  genes.info_df<-as_tibble(data.frame(chr,start, end, gene_name, data))
  
  #remove those probes without annotated genes
  genes.info_df<-genes.info_df[genes.info_df$gene_name != ".", ]
  
  
  #operations to drop the first three chracters in the chromosome column
  genes.info_df$chr<-unfactor(genes.info_df$chr)
  
  dropchr <- function(x) {  
    substr(x,4,4)
    }
  
  genes.info_df$chr<-sapply(genes.info_df$chr, dropchr)
  
  
  
  ########################################################
  #4. Finding overlaps
  ########################################################
  
  cat("\nFinding overlaps between the annotation and our methylation probes.\n")
  
  #make GRanges from annotation
  colnames(genes) <- c("chr","start", "end", "GeneSymbol")
  genes$start<-genes$start-3000
  genes_GR <- makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
  
  #make GRanges object from the tibble
  df_GR <- makeGRangesFromDataFrame(genes.info_df[,c(1,2,3)])#we only take chr, start and end to do the overlap
  
  #overlap the annotation with the cpgs
  hits <- findOverlaps(df_GR,genes_GR, type="within") #hits found
  
  cat("\nCreating data frame of useful and annotated cpgs.\n")
  
  df_ann <- cbind(genes.info_df[queryHits(hits),],genes[subjectHits(hits),])
  drops <- c("chr","start","end","gene_name","position_to_TSS")
  clean_df<-df_ann[ , !(names(df_ann) %in% drops)]
  clean_df<-clean_df[c(24,1:23)]
  
  ########################################################
  #4. Agregation of beta value means
  ########################################################
  
  #remove the NA that are found across all columns
  clean_df<-clean_df[complete.cases(clean_df), ]
  
  cat("\nAggregating cpg values by gene.\n")
  
  agg_df <- aggregate(x = clean_df[2:ncol(clean_df)], 
                                 by = list(clean_df$GeneSymbol), 
                                 FUN = mean)
  
  rownames(agg_df)<-agg_df$Group.1
  drops <- "Group.1"
  agg_df<-agg_df[ , !(names(agg_df) %in% drops)]
  
  totalmean<-apply(agg_df, 1, mean)
  agg_df<-cbind(totalmean, agg_df)
  
  return(agg_df)
}


########################################################
#Testing the function
########################################################
LUAD.pts <- task1("TCGA-LUAD")
KIRC.pts <- task1("TCGA-KIRC")
HNSC.pts <- task1("TCGA-HNSC")
STAD.pts <- task1("TCGA-STAD")
LUSC.pts <- task1("TCGA-LUSC")
KICH.pts <- task1("TCGA-KICH")
SKCM.pts <- task1("TCGA-SKCM")
KIRP.pts <- task1("TCGA-KIRP")
ESCA.pts <- task1("TCGA-ESCA")

LUAD.meth <- task4("TCGA-LUAD",  LUAD.pts)
KIRK.meth <- task4("TCGA-KIRC",  KIRC.pts)
HNSC.meth <- task4("TCGA-HNSC",  HNSC.pts)
STAD.meth <- task4("TCGA-STAD",  STAD.pts )
LUSC.meth <- task4("TCGA-LUSC",  LUSC.pts)
KICH.meth <- task4("TCGA-KICH",  KICH.pts)
SKCM.meth <- task4("TCGA-SKCM",  SKCM.pts )
KIRP.meth <- task4("TCGA-KIRP",  KIRP.pts)
ESCA.meth <- task4("TCGA-ESCA",  ESCA.pts)