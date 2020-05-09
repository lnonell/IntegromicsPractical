#Task 8 Oscar Paniello Pérez-Hickman
#Purpose: Plot results of MFA, mixOmics, together with rawdata (mRNA, CN and methylation data) and significative correlation data.
#Filter previously each raw data set to get the 10% of the genes with the highest sd. The plot should be informative.
#input: data.frames obtained as output files in tasks 1 to 5. 100 most correlated variables obtained in tasks 6 and 7. Path to store Circos plot
#output: Circos plot in specified format. 



task8<-function(df_samples, mrna, cnv, meth, cor.cnv, cor.cnv.sig = 3, cor.meth, cor.meth.sig = 3, sd_mrna = 0.01, sd_cnv = 1, sd_meth = 0.01,path = getwd()){
  
  #library
  suppressPackageStartupMessages(library(OmicCircos))
  suppressPackageStartupMessages(library(biomaRt))
  suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(compare))
  suppressPackageStartupMessages(library(TCGAbiolinks))
  suppressPackageStartupMessages(library(edgeR))
  suppressPackageStartupMessages(library(SummarizedExperiment))
  
  # filter 10% genes by highest sd
  cat("\nfilter 10% genes by highest sd.\n")
  
  sd <- apply(mrna,1,sd)
  sd10 <- head(sort(sd,decreasing = TRUE), round(nrow(mrna)*sd_mrna))
  exp_final <- mrna[names(sd10),]
  
  sd <- apply(cnv,1,sd)
  sd10 <- head(sort(sd,decreasing = TRUE), round(nrow(cnv)*sd_cnv))
  cnv_final <- cnv[names(sd10),]
  
  sd <- apply(meth,1,sd)
  sd10 <- head(sort(sd,decreasing = TRUE), round(nrow(meth)*sd_meth))
  meth_final <- meth[names(sd10),]
  
  
  cat("\nDownloading annotation of human genes.\n")
  
  #Set up an gene annotation template to use
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genes <- genes(txdb)
  ensembl <- useMart("ensembl")
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
  
  #obtain chr, start, and HGNC name of all the genes annotated in hg38
  annot_df <- getBM(attributes = c("chromosome_name","start_position","end_position","hgnc_symbol"), filters = "entrezgene_id", values = genes$gene_id, mart = mart)
  
  annot_df<- annot_df[annot_df[,4]!="" & annot_df[,1] %in% c(1:22,"X","Y"),] #remove those not annotated 
  annot_df <- annot_df[order(annot_df[,2]),] #order by start position
  genes <- annot_df[order(annot_df[,1]),] #order by chromosome
  
  genes<-genes %>% distinct(hgnc_symbol, .keep_all = TRUE) #we select just one copy of the gene with it starting and end positon
  rownames(genes)<-genes$hgnc_symbol
  
  cat("\nChecking the common genes between the df and the annotation.\n")
  #Obtain just the positions of the wanted genes in common with the annotation gene list.
  #Exp
  common.gene.exp <- intersect(row.names(genes), row.names(exp_final))
  gene.pos.exp <- genes[common.gene.exp,]
  drops <- c("hgnc_symbol")
  gene.pos.exp <- gene.pos.exp[ , !(names(gene.pos.exp) %in% drops)]
  exp_final <- exp_final[common.gene.exp,]
  
  #CNV
  common.gene.cnv <- intersect(row.names(genes), row.names(cnv_final))
  gene.pos.cnv <- genes[common.gene.cnv,]
  drops <- c("hgnc_symbol")
  gene.pos.cnv <- gene.pos.cnv[ , !(names(gene.pos.cnv) %in% drops)]
  cnv_final <- cnv_final[common.gene.cnv,]
  
  cnv.m <- merge(gene.pos.cnv, cnv_final, by=0, all=TRUE)  # merge by row names (by=0 or by="row.names")
  cnv.m[is.na(cnv.m)] <- 0                 # replace NA values
  cnv <- cnv.m[,-1]
  rownames(cnv) <- cnv.m[,1]
  
  #Meth
  common.gene.meth <- intersect(row.names(genes), row.names(meth_final))
  gene.pos.meth <- genes[common.gene.meth,]
  drops <- c("hgnc_symbol")
  gene.pos.meth <- gene.pos.meth[ , !(names(gene.pos.meth) %in% drops)]
  meth_final <- meth_final[common.gene.meth,]
  
  #Correlation exp-cnv
  common.gene.corr.cnv <- intersect(row.names(genes), row.names(cor.cnv))
  gene.pos.corr.cnv <- genes[common.gene.corr.cnv,]
  drops <- c("hgnc_symbol")
  gene.pos.corr.cnv <- gene.pos.corr.cnv[ , !(names(gene.pos.corr.cnv) %in% drops)]
  
  #-log10(pvalue) to make a better visual correlation in the circos
  logpvalue.cnv.corr <- -1 * log10(cor.cnv[,3]) 
  cor.cnv[,3] <- logpvalue.cnv.corr
  
  #prepare data for the circosplot
  pvalue.exp.cnv <- merge(gene.pos.corr.cnv, cor.cnv, by="row.names",all.x=TRUE)
  rownames(pvalue.exp.cnv)<-pvalue.exp.cnv[,1]
  pvalue.exp.cnv <- pvalue.exp.cnv[,c(-1,-5,-6)]
  
  mask.exp.cnv <- pvalue.exp.cnv[,"pval.adj"]>cor.cnv.sig
  sig.pvalue.exp.cnv <- pvalue.exp.cnv[mask.exp.cnv,]
  sig.pvalue.exp.cnv[,5] <- rownames(sig.pvalue.exp.cnv)
  
  #Correlation exp-meth
  common.gene.corr.meth <- intersect(row.names(genes), row.names(cor.meth))
  gene.pos.corr.meth <- genes[common.gene.corr.meth,]
  drops <- c("hgnc_symbol")
  gene.pos.corr.meth <- gene.pos.corr.meth[ , !(names(gene.pos.corr.meth) %in% drops)]
  
  #-log10(pvalue) to make a better visual correlation in the circos
  logpvalue.meth.corr <- -1 * log10(cor.meth[,3])
  cor.meth[,3] <- logpvalue.meth.corr
  
  #prepare data for the circosplot
  pvalue.exp.meth <- merge(gene.pos.corr.meth, cor.meth, by="row.names",all.x=TRUE)
  rownames(pvalue.exp.meth)<-pvalue.exp.meth[,1]
  pvalue.exp.meth <- pvalue.exp.meth[,c(-1,-5,-6)]
  mask.exp.meth <- pvalue.exp.meth[,"pval.adj"]>cor.meth.sig
  sig.pvalue.exp.meth <- pvalue.exp.meth[mask.exp.meth,]
  sig.pvalue.exp.meth[,5] <- rownames(sig.pvalue.exp.meth)
  
  circ.exp <- data.frame("chr"=gene.pos.exp$chromosome_name,"Start"=as.integer(gene.pos.exp$start_position),
                         "End"=as.integer(gene.pos.exp$end_position),log2(exp_final+1),row.names=NULL)
  circ.meth <- data.frame("chr"=gene.pos.meth$chromosome_name,"Start"=as.integer(gene.pos.meth$start_position),
                          "End"=as.integer(gene.pos.meth$end_position),meth_final,row.names=NULL)
  
  cat("\nPlotting the circos as pdf in the working directory.\n")
  pdf("circosplot.pdf")
  
  options(stringsAsFactors = FALSE)
  par(mar=c(2, 2, 2, 2));
  plot(c(1,800), c(1,800), type="n", axes=F, xlab="", ylab="", main="");
  circos(R=400, cir="hg19", W=4,   type="chr", print.chr.lab=T, scale=T);
  circos(R=300, cir="hg19", W=100,  mapping=circ.exp,   col.v=4,    type="heatmap2",B=FALSE, cluster=FALSE, col.bar=F, lwd=0.1, col="blue");
  circos(R=120, cir="hg19", W=80,  mapping=cnv,   col.v=4,   type="ml3", B=FALSE, lwd=1, cutoff=0);
  circos(R=200, cir="hg19", W=100,  mapping=circ.meth,   col.v=4, type="heatmap2",B=FALSE, cluster=FALSE, col.bar=F, lwd=0.1, col="blue")
  colors <- rainbow( 1 , start=0, end = 1, alpha=0.5 )
  circos(R=90, cir="hg18", W=40,  mapping=pvalue.exp.cnv,  col.v=4,    type="s",   B=FALSE, lwd=0.1, col=colors, cex=0.2)
  colors <- rainbow( 1 , start = 0.7, end = 1, alpha=0.5 )
  circos(R=90, cir="hg18", W=40,  mapping=pvalue.exp.meth,  col.v=4,    type="s",   B=FALSE, lwd=0.1, col=colors, cex=0.2)
  if (nrow(sig.pvalue.exp.cnv) > 0) {
    circos(R=120, cir="hg18", W=40,  mapping=sig.pvalue.exp.cnv,  col.v=5,    type="label2",   B=FALSE, lwd=1, col='red', cex=0.2)
  }
  if (nrow(sig.pvalue.exp.meth) > 0) {
    circos(R=120, cir="hg18", W=40,  mapping=sig.pvalue.exp.meth,  col.v=5,    type="label2",   B=FALSE, lwd=1, col='blue', cex=0.2)
  }
  
  dev.off() 
  
  
  return()
}

##############
#Testing
##############
task8(LUAD.pts)
task8(KIRC.pts)
task8(HNSC.pts)
task8(STAD.pts)
task8(LUSC.pts)
task8(KICH.pts)
task8(SKCM.pts, mrna = SKCM.exp, cnv = SKCM.cnv, meth = SKCM.meth, cor.cnv = SKCM.cnv.corr, cor.meth = SKCM.meth.corr, path = getwd())
task8(df_samples = KIRP.pts,mrna = KIRP.exp,cnv = KIRP.cnv,meth = KIRP.meth,cor.cnv = KIRP.cnv.corr,cor.meth = KIRP.meth.corr)
task8(ESCA.pts)

