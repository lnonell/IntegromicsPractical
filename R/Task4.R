#Task 4: Cristina Tuni i Dominguez
#Purpose: Prepare methylation data to be downloaded and further treated.
#input: Data frame with the prepared clinical data from task 1, plus the cancer code from TCGA.
#output: Data frame with samples in columns and genes (HUGO Symbol) in rows containing the beta values mean of each gene and sample.

task4<-function(cancer, df_samples, total.mean=FALSE){
  
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

  
  # get expression matrix
  data<-assay(TCGA.meth)
  colnames(data)<-substr(colnames(data),1,12)
  
  #removal of duplicates (will remove the first duplicated column)
  data <- data[, !duplicated(colnames(data))]
  
  #get genes information
  genes.info<- rowRanges(TCGA.meth)
  
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
    substring(x,4)
    }
  
  genes.info_df$chr<-sapply(genes.info_df$chr, dropchr)
  
  
  
  ########################################################
  #4. Finding overlaps
  ########################################################
  
  cat("\nFinding overlaps between the annotation and our methylation probes.\n")
  
  #make GRanges from annotation
  genes$end_position<-genes$start_position
  genes$start_position<-genes$start_position-3000
  colnames(genes) <- c("chr","start", "end", "GeneSymbol")
  genes_GR <- makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
  
  #make GRanges object from the tibble
  df_GR <- makeGRangesFromDataFrame(genes.info_df[,c(1,2,3)])#we only take chr, start and end to do the overlap
  
  #overlap the annotation with the cpgs
  hits <- findOverlaps(df_GR,genes_GR, type="within") #hits found
  
  cat("\nCreating data frame of useful and annotated cpgs.\n")
  
  df_ann <- cbind(genes.info_df[queryHits(hits),],genes[subjectHits(hits),])
  drops <- c("chr","start","end","gene_name","position_to_TSS")
  clean_df<-df_ann[ , !(names(df_ann) %in% drops)]
  
  
  
  #Function to move the last column to the first one, by StackOverflow user: A5C1D2H2I1M1N2O1R2T1
  moveme <- function (invec, movecommand) {
    movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                   ",|\\s+"), function(x) x[x != ""])
    movelist <- lapply(movecommand, function(x) {
      Where <- x[which(x %in% c("before", "after", "first", 
                                "last")):length(x)]
      ToMove <- setdiff(x, Where)
      list(ToMove, Where)
    })
    myVec <- invec
    for (i in seq_along(movelist)) {
      temp <- setdiff(myVec, movelist[[i]][[1]])
      A <- movelist[[i]][[2]][1]
      if (A %in% c("before", "after")) {
        ba <- movelist[[i]][[2]][2]
        if (A == "before") {
          after <- match(ba, temp) - 1
        }
        else if (A == "after") {
          after <- match(ba, temp)
        }
      }
      else if (A == "first") {
        after <- 0
      }
      else if (A == "last") {
        after <- length(myVec)
      }
      myVec <- append(temp, values = movelist[[i]][[1]], after = after)
    }
    myVec
  }
  
  
  clean_df<-clean_df[moveme(names(clean_df), "GeneSymbol first")]
  
  ########################################################
  #4. Agregation of beta value means
  ########################################################
  
  cat("\nAggregating cpg values by gene.\n")
  
  mean.na<-function(x){
    mean(x, na.rm = TRUE)
  }
  
  agg_df <- aggregate(x = clean_df[2:ncol(clean_df)], 
                                 by = list(clean_df$GeneSymbol), 
                                 FUN = mean.na)
  
  agg_df<-na.omit(agg_df)
  rownames(agg_df)<-agg_df$Group.1
  drops <- "Group.1"
  agg_df<-agg_df[ , !(names(agg_df) %in% drops)]
  colnames(agg_df)<-gsub("\\.","-", colnames(agg_df))
  
  if (total.mean==FALSE){
    return(agg_df)
  } else {
    totalmean<-apply(agg_df, 1, mean)
    agg_df<-cbind(totalmean, agg_df)
    return(agg_df)
  }
  
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
KIRC.meth <- task4("TCGA-KIRC",  KIRC.pts)
HNSC.meth <- task4("TCGA-HNSC",  HNSC.pts)
STAD.meth <- task4("TCGA-STAD",  STAD.pts )
LUSC.meth <- task4("TCGA-LUSC",  LUSC.pts)
KICH.meth <- task4("TCGA-KICH",  KICH.pts)
SKCM.meth <- task4("TCGA-SKCM",  SKCM.pts )
KIRP.meth <- task4("TCGA-KIRP",  KIRP.pts)
ESCA.meth <- task4("TCGA-ESCA",  ESCA.pts)


cancer<-"TCGA-KIRC"
df_samples<-KIRC.pts
