#Task 7 Gerard Temprano
#Purpose
#input
#output
task7<-function(){# Task 7: Gerard Temprano Sagrera
  # Purpose: Apply mixOmics to the data set comparing stage i vs stage iv.
  # Input: Data-frames with mRNA, CN and methylation data and destination path. 
  # Output: Inoformative graphs and the 100 most correlated variables.
  
  task7<-function(Cancer_table, mRNA, CNV, Methylation, Destination_path){
    
    if (class(mRNA) != "data.frame") {
      stop(print("mRNA is not a data.frame"))}
    
    if (class(CNV) != "data.frame") {
      stop(print("CNV is not a data.frame"))}
    
    if (class(Methylation) != "data.frame") {
      stop(print("Methylation is not a data.frame"))}
  }
  
  task7(LUSC.table, mRNA, breast.TCGA, mRN)
  
  
  # Select the 10% of the genes with the highest standard deviation for each of the 3 data types:
  mRNA
  SD<-apply(mRNA,1,sd)
  top.10sd<-head(sort(SD,decreasing=TRUE),round(length(mRNA)/10))
  
  mRNA.f <- mRNA[names(top.10sd),]
  
  
  CNV
  SD<-apply(CNV,1,sd)
  top.10sd<-head(sort(SD,decreasing=TRUE),round(length(CNV)/10))
  
  CNV.f <- CNV[names(top.10sd),]
  
  
  Methylation
  SD<-apply(Methylation,1,sd)
  top.10sd<-head(sort(SD,decreasing=TRUE),round(length(Methylation)/10))
  
  Methylation.f <- Methylation[names(top.10sd),]
  
  
  # Compare stage i vs stage iv:
  
  # Cancer_table
  stage_iv <- grep("stage iv", LUSC.table$tumor_stage)
  stage_i <- length(grep("stage i", LUSC.table$tumor_stage)) - stage_iv
  
  cond <- c(rep("stage i", length(stage_i)), rep("stage iv", length(stage_iv)))
  Y=as.factor(cond)
  
  
  data = list(mRNA = mRNA.f, 
              CNA = CNV.f,
              meta=Methylation.f)
  
  design = matrix(0.1, ncol = length(data), nrow = length(data), 
                  dimnames = list(names(data), names(data)))
  
  # To avoid analyzing an experiment against itself
  diag(design) = 0
  
  design
  
  
  # Classify a discrete outcome through the integration of multiple datasets selecting features from each data set, using PLS-DA model:
  sgccda.res = block.splsda(X = data, Y = Y, ncomp = 2, design = design)
  
  set.seed(123)
  # Evaluate performance of the PLS-DA model:
  perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)
  
  pdf("Diablo Performance")
  plot(perf.diablo) # Components that exaplain most of the variability 
  dev.off()
  
  perf.diablo$choice.ncomp$WeightedVote
  ncomp = 2
  
  # Tune the keepX parameters obtained with the block.splsda function:
  set.seed(1234) 
  s.keep <- sort(sample(1:150,15))
  test.keepX = list (rna = s.keep, dna =s.keep)
  
  # Computes M-fold scores based on a user-input grid to determine the optimal parsity parameters 
  tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp,
                                test.keepX = test.keepX, design = design,
                                validation = 'Mfold', folds = 10, nrepeat = 1,
                                cpus = 4, dist = "centroids.dist")
  
  list.keepX = tune.TCGA$choice.keepX
  tune.TCGA$choice.keepX
  
  list.keepX
  
  ###### Apply method with tunned params
  sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                            keepX = list.keepX, design = design)
  
  sgccda.res$design
  
  ###### Results
  head(selectVar(sgccda.res)$dna$value,10)
  head(selectVar(sgccda.res)$rna$value,10)
  
  ###### Plots
  pdf("DIABLO.Pollack.pdf")
  plotDiablo(sgccda.res, ncomp = 1)
  
  plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
  plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
          pch = c(16, 17), cex = c(2,2), col = c('darkorchid', 'lightgreen'))
  
  plotLoadings(sgccda.res, comp = 2, contrib = 'max', method = 'median')
  
  cimDiablo(sgccda.res)
  
  dev.off()
  
  
  
  
  
  
  
  
  
  
  return()
}