cancer = "TCGA-LUAD"
df_samples = LUAD.pts
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
  
  ########################################################
  #1. Downloading methylation data
  ########################################################
  
  query.meth <- GDCquery(project = cancer, 
                         data.category = "DNA Methylation",
                         platform = "Illumina Human Methylation 450", 
                         barcode = df_samples$barcode)
  GDCdownload(query.meth)
  TCGA.meth<- GDCprepare(query=query.meth, 
                         save=T, 
                         save.filename = "DNAmeth.rda", 
                         summarizedExperiment = T)
  
  # get expression matrix
  data<-assay(TCGA.meth)
  
  #get genes information
  genes.info<- rowRanges(TCGA.meth)
  
  #get sample information
  sample.info<- colData(TCGA.meth)
  
  ########################################################
  #2. Annotation
  ########################################################
  
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
  
  #check if the data object and the genes.info object have the same order (they should)
  all(genes.info@ranges@NAMES==rownames(data)) #TRUE
  
  #from the previous GRanges, obtain useful dataframe
  
  chr<-genes.info@seqnames
  
  start<-as.double(genes.info@ranges@start)
  end<-as.double(start_position+1)
  
  gene_name<-genes.info@elementMetadata@listData$Gene_Symbol
  i=0
  for (i in 1:length(gene_name)){
    gene_name[i]<-strsplit(gene_name[i], "[;]")[[1]][1]
  }
  
  position_to_TSS<-genes.info@elementMetadata@listData$Position_to_TSS
  i=0
  for (i in 1:length(position_to_TSS)){
    position_to_TSS[i]<-strsplit(position_to_TSS[i], "[;]")[[1]][1]
  }
  
  
  genes.info_df<-as_tibble(data.frame(chr,start, end, gene_name, position_to_TSS, data))
  
  #remove those probes without annotated genes
  genes.info_df<-genes.info_df[genes.info_df$gene_name != ".", ]
  
  
  #operations to drop the first three chracters in the chromosome column
  genes.info_df$chr<-unfactor(genes.info_df$chr)
  
  dropchr <- function(x) {  
    substr(x,4,4)}
  
  genes.info_df$chr<-sapply(genes.info_df$chr, dropchr)
  
  #transform the positions to numeric type
  genes.info_df$position_to_TSS<-as.numeric(as.character(genes.info_df$position_to_TSS))
  
  
  ########################################################
  #4. Finding overlaps
  ########################################################
  
  #make GRanges from annotation
  colnames(genes) <- c("chr","start", "end", "GeneSymbol")
  genes$start<-genes$start-3000
  genes_GR <- makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
  
  #make GRanges object from the tibble
  df_GR <- makeGRangesFromDataFrame(genes.info_df[,c(1,2,3)])#we only take chr, start and end to do the overlap
  
  #overlap the annotation with the cpgs
  hits <- findOverlaps(df_GR,genes_GR, type="within") #hits found
  
  df_ann <- cbind(genes.info_df[subjectHits(hits),],genes[queryHits(hits),]) 
  
  
  #aggregated data frame
  first_sample<-aggregate(TCGA.05.4390.01A.02D.1756.05 ~ gene_name, FUN = mean, data=genes.info_df)
  
  #remove the NA that are found across all columns
  genes.info_df<-genes.info_df[complete.cases(genes.info_df), ]
  
  #remove those probes without annotated genes
  genes.info_df<-genes.info_df[genes.info_df$gene_name != ".", ]
  
  #transform the positions to numeric type
  genes.info_df$position_to_TSS<-as.numeric(as.character(genes.info_df$position_to_TSS))
  
  
  
  
  return()
}
