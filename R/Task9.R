#Task 9
#Purpose Functional analysis using KEGG
#input: correlated variables obtained in tasks 6 and 7. Path to store files.
#output top KEGG path painted with pathview. Barplot with most significant gene sets.

#### Important message on imputs
#Go to task6, final line (below return & before the keys) and add the following code:
#MFA6 <- data.frame(PC1, PC2)
#Rerun TASK6 with your cancer type.

#Same for  task 7, final line (below return & before the keys) and add the following code:
#DIABLO7 <- list(data.table1, data.table2)
#Rerun TASK7 with your cancer type.
#Proceed when objects MFA6 & DIABLO7 are created.

task9<-function(MFA6, DIABLO7){
  
  #####
  # Load libraries
  suppressPackageStartupMessages(library(GSEABase))
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  suppressPackageStartupMessages(library(KEGG.db))
  suppressPackageStartupMessages(library(GOstats))
  suppressPackageStartupMessages(library(pathview))
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(enrichplot))
  
  #############################################################################################
  #### Functional analysis of MFA results (Task6)
  ############################################################################################
  
 
  # Go to task 6, final line (below return) and add the following code:
  cat("Go to task 6, final line (below return & before the keys) and add the following code:\n MFA6 <- data.frame(PC1, PC2) \n rerun TASK6 with your cancer type. \n Proceed when object MFA6 is created.")
  stopifnot(is.data.frame(MFA6))
  cat("Go to task 7, final line (below return & before the keys) and add the following code:\n DIABLO7 <- list(data.table1, data.table2) \n rerun TASK7 with your cancer type. \n Proceed when object DIABLO7 is created.")
  stopifnot(is.list(DIABLO7))
  

  # Remove version from gene IDs
  MFA6_PC1 <- MFA6$PC1
  MFA6_PC1_nv <- gsub('.{2}$', '', MFA6_PC1)
  class(MFA6_PC1_nv)
  MFA6_PC2 <- MFA6$PC2
  MFA6_PC2_nv <- gsub('.{2}$', '', MFA6_PC2)
  
  #Unify the most correlated genes from obtained using MFA
  MFA6_nv <- c(MFA6_PC1_nv, MFA6_PC2_nv)
  
  # Transform gene IDs to Entrez IDs
  MFA6_nv_entrez <- mapIds(org.Hs.eg.db, MFA6_nv, 'ENTREZID', 'SYMBOL')
  
  #Create KEGG GeneSetCollection based in symbol
  frame = toTable(org.Hs.egPATH)
  frame$symbol<-getSYMBOL(frame$gene_id,data = 'org.Hs.eg')
  keggframeData = data.frame(frame$path_id, frame$gene_id)
  keggFrame=KEGGFrame(keggframeData,organism="Homo sapiens")
  gsc.KEGG <- GeneSetCollection(keggFrame, setType = KEGGCollection())
  
  # Establish the parameters for the KEGG pathway enrichment analysis
  cutoff=0.01
  KEGG.params <- GSEAKEGGHyperGParams(name="KEGG",  geneSetCollection=gsc.KEGG, geneIds = MFA6_nv_entrez, universeGeneIds = NULL,  pvalueCutoff = cutoff, testDirection = "over")
  KEGG.results<-hyperGTest(KEGG.params)
  KEGG.results
  summary(KEGG.results)
  
  # Create results repository
  dir.create(file.path(getwd(), "./kegg_results"), showWarnings = FALSE)
  dir.create(file.path(getwd(), "./kegg_results/MFA"), showWarnings = FALSE)
  
  #Export the pathview representation
  path.id<-summary(KEGG.results)$KEGGID
  pv.out <- pathview(gene.data=MFA6_nv_entrez, gene.idtype="entrez",pathway.id=path.id, species="hsa", kegg.dir = "./kegg_results/", out.suffix="MFA.kegg")
  str(pv.out)

  #Pathway dictionary. To facilitate KEGGID interpretation
  dict1 <- data.frame(KEGGID = summary(KEGG.results)$KEGGID, Pathway = summary(KEGG.results)$Term)
  write.table(dict1, "./kegg_results/MFA/keggid2pathway_dictionary1.tsv", quote = F, sep = "\t")
  
  #KEGG enrichment results (Hypergeometric test)
  results1 <- data.frame(summary(KEGG.results))
  write.table(results1, "./kegg_results/MFA/kegg_summary_MFA.tsv", quote = F, sep = "\t")
  
  ## 
  # Perform enrichment analysis using clusterProfiler
  ans.kegg <- enrichKEGG(gene = MFA6_nv_entrez,
                         organism = 'hsa',
                         pvalueCutoff = 0.1)
  
  # Export the barplot with the most significant gene sets
  pdf("./kegg_results/MFA/barplot_MFA_enriched_genes.pdf")
  barplot(ans.kegg, showCategory=20, title = "Enriched Gene Sets")
  dev.off()
  
  #############################################################################################
  #### Functional analysis of DIABLO results (Task7)
  ############################################################################################

  # Remove version from gene IDs
  DIABLO7[[1]]$Genes <- as.vector(DIABLO7[[1]]$Genes)
  DIABLO7[[2]]$Genes <- as.vector(DIABLO7[[2]]$Genes)
  DIABLO7_genes <- c(DIABLO7[[1]]$Genes, DIABLO7[[2]]$Genes)
  DIABLO7_genes_nv <- gsub('.{2}$', '', DIABLO7_genes)
  
  # Transform gene IDs to Entrez IDs
  DIABLO7_nv_entrez <- mapIds(org.Hs.eg.db, DIABLO7_genes_nv, 'ENTREZID', 'SYMBOL')
  
  #Create KEGG GeneSetCollection based in symbol
  frame = toTable(org.Hs.egPATH)
  frame$symbol<-getSYMBOL(frame$gene_id,data = 'org.Hs.eg')
  keggframeData = data.frame(frame$path_id, frame$gene_id)
  keggFrame=KEGGFrame(keggframeData,organism="Homo sapiens")
  gsc.KEGG <- GeneSetCollection(keggFrame, setType = KEGGCollection())
  
  # Establish the parameters for the KEGG pathway enrichment analysis
  cutoff=0.01
  KEGG.params <- GSEAKEGGHyperGParams(name="KEGG",  geneSetCollection=gsc.KEGG, geneIds = DIABLO7_nv_entrez, universeGeneIds = NULL,  pvalueCutoff = cutoff, testDirection = "over")
  KEGG.results<-hyperGTest(KEGG.params)
  KEGG.results
  summary(KEGG.results)
  
  # Create results repository
  dir.create(file.path(getwd(), "./kegg_results/DIABLO"), showWarnings = FALSE)
  
  #Export the pathview representation
  path.id<-summary(KEGG.results)$KEGGID
  pv.out <- pathview(gene.data=DIABLO7_nv_entrez, gene.idtype="entrez",pathway.id=path.id, species="hsa", kegg.dir = "./kegg_results/DIABLO", out.suffix="DIABLO.kegg")
  str(pv.out)
  
  #Pathway dictionary. To facilitate KEGGID interpretation
  dict2 <- data.frame(KEGGID = summary(KEGG.results)$KEGGID, Pathway = summary(KEGG.results)$Term)
  write.table(dict2, "./kegg_results/DIABLO/keggid2pathway_dictionary2.tsv", quote = F, sep = "\t")
  
  #KEGG enrichment results
  results2 <- data.frame(summary(KEGG.results))
  write.table(results2, "./kegg_results/DIABLO/kegg_summary_DIABLO.tsv", quote = F, sep = "\t")

  ## 
  # Perform enrichment analysis using clusterProfiler
  ans.kegg <- enrichKEGG(gene = DIABLO7_nv_entrez,
                         organism = 'hsa',
                         pvalueCutoff = 0.1)
  
  # Export the barplot with the most significant gene sets
  pdf("./kegg_results/DIABLO/barplot_DIABLO_enriched_genes.pdf")
  barplot(ans.kegg, showCategory=20, title = "Enriched Gene Sets")
  dev.off()
  
  
  
  return(cat("Find the results at the folder kegg_results"))  
}
