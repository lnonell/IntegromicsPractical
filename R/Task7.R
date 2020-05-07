getwd()
# Task 7: Gerard Temprano Sagrera
# Purpose: Apply mixOmics to the data set comparing stage i vs stage iv.
# Input: Data-frames with mRNA, CN and methylation data and destination path. 
# Output: Inoformative graphs and the 100 most correlated variables.

task7<-function(Cancer_table, mRNA, CNV, Methylation){
  
  ################
  #Load libraries#
  ################
  
  suppressPackageStartupMessages(library(mixOmics))
  
  
  if(missing(Methylation)){
    print("Methylation data not introduced")
    
    #####################
    #STEP 1: PREPARATION#
    #####################
    
    # Check if the objects are data frames #
    
    if (class(Cancer_table) != "data.frame") {
      stop(print("Cancer_table is not a data.frame"))}
    
    if (class(mRNA) != "data.frame") {
      stop(print("mRNA is not a data.frame"))}
    
    if (class(CNV) != "data.frame") {
      stop(print("CNV is not a data.frame"))}
    
    
    # Set different gene names for each datatable #
    
    rownames(mRNA) <- paste(rownames(mRNA),"r",sep=".")
    rownames(CNV) <- paste(rownames(CNV),"c",sep=".") 
    
    
    ## Set samples names prepared for mixOmics with the same format and order ##
    
    
    # Change separations in samples names #
    
    colnames(mRNA) <- gsub("\\.", "-",colnames(mRNA))
    colnames(CNV) <- gsub("\\.","-",colnames(CNV))
    
    
    # Order to allow stage comparisons #
    
    mRNA <- mRNA[, order(match(colnames(mRNA), Cancer_table[,1]))]
    CNV <- CNV[, order(match(colnames(CNV), Cancer_table[,1]))]
    
    
    if (all(sapply(list(colnames(mRNA),colnames(CNV), Cancer_table[,1]), function(x) x == Cancer_table[,1])) == FALSE){
      stop(print("Samples are not the same"))}
    
    
    # Set inputs as numeric data frames #
    
    if (class(mRNA[1,1]) != c("numeric")){
      mRNA2 <- as.data.frame(sapply(mRNA, as.numeric))
      rownames(mRNA2) <- rownames(mRNA)   
      mRNA <- mRNA2
      print("Expression dataset changed to numeric data frame")}
    
    if (class(CNV[1,1]) != c("numeric")){
      CNV2 <- as.data.frame(sapply(CNV, as.numeric))
      rownames(CNV2) <- rownames(CNV)
      CNV <- CNV2
      print("CNV dataset changed to numeric data frame")}    
    
    
    # Select the 10% of the genes with the highest standard deviation for each of the 3 data types #
    
    SD<-apply(mRNA,1,sd)
    top.10sd<-head(sort(SD,decreasing=TRUE),round(length(SD)/10))
    mRNA.f <- mRNA[names(top.10sd),]
    
    SD<-apply(CNV,1,sd)
    top.10sd<-head(sort(SD,decreasing=TRUE),round(length(SD)/10))
    CNV.f <- CNV[names(top.10sd),]
    
    
    #######################
    #STEP 2: APPLY  DIABLO#
    #######################
    
    
    
    # List of data sets #
    
    data = list(mRNA = t(mRNA.f),
                CNV = t(CNV.f))
    
    
    # Set Y to compare stage i vs stage iv #
    
    Y=as.factor(Cancer_table$tumor_stage)
    
    design = matrix(0.1, ncol = length(data), nrow = length(data), 
                    dimnames = list(names(data), names(data)))
    lapply(data,dim)
    Y
    diag(design) = 0 # To avoid analyzing an experiment against itself
    design
    
    
    sgccda.res = block.splsda(X = data, Y = Y, ncomp = 5, design = design) # Classifies a discrete outcome through the integration of 
    # multiple datasets selecting features from each data set, using a PLS-DA model
    
    
    # Evaluate performance of the PLS-DA model #
    
    set.seed(12345)
    
    perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)
    
    print("Saving results in Diablo Results folder in your working directory")
    dir.create(file.path(getwd(), "./Diablo Results"), showWarnings = FALSE)
    png("./Diablo Results/Diablo first model performance.png")
    plot(perf.diablo) # Components that exaplain most of the variability
    dev.off()
    
    ncomp = 2
    
    # Tune the keepX parameters obtained with the block.splsda function #
    
    s.keep <- sort(sample(1:length(rownames(Cancer_table)),5))
    test.keepX = list(mRNA = s.keep, CNV =s.keep)
    
    
    # Computes M-fold scores based on a user-input grid to determine the optimal parsity parameters #
    
    tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = 2,
                                  test.keepX = test.keepX, design = design,
                                  validation = 'Mfold', folds = 10, nrepeat = 1,
                                  cpus = 4, dist = "centroids.dist")
    
    list.keepX = tune.TCGA$choice.keepX
    
    
    # Apply method with tunned params #
    sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                              keepX = list.keepX, design = design)
    
    sgccda.res$design
    
    
    ################
    #STEP3: RESULTS#
    ################
    
    # Component 1 #
    
    varm1 <- head(selectVar(sgccda.res, comp = 1)$mRNA$value,100)
    varcnv1 <- head(selectVar(sgccda.res, comp = 1)$CNV$value,100)
    
    data.table1 <- rbind(varm1, varcnv1)
    data.table1 <- cbind(data.table1, rownames(data.table1))
    data.table1 <- data.table1[order(abs(data.table1$value.var), decreasing = TRUE),]
    colnames(data.table1) <- c("value.var1", "Genes")
    
    
    # Component 2 #
    
    varm2 <- head(selectVar(sgccda.res, comp = 2)$mRNA$value,100)
    varcnv2 <- head(selectVar(sgccda.res, comp = 2)$CNV$value,100)
    
    data.table2 <- rbind(varm2, varcnv2)
    data.table2 <- cbind(data.table2, rownames(data.table2))
    data.table2 <- data.table2[order(abs(data.table2$value.var), decreasing = TRUE),]
    colnames(data.table2) <- c("value.var2", "Genes")
    
    
    if (nrow(data.table1) > 100) {
      data.table1 <- data.table1[1:100,]}
    
    if (nrow(data.table2) > 100) {
      data.table2 <- data.table2[1:100,]}
    
    ###############
    #STEP4: PLOTS#
    ###############
    
    png("./Diablo Results/Diablo correlation first component.png")
    plotDiablo(sgccda.res, ncomp = 1)
    dev.off()
    
    png("./Diablo Results/Diablo correlation second component.png")
    plotDiablo(sgccda.res, ncomp = 2)
    dev.off()
    
    png("./Diablo Results/Arrow plot samples.png")
    plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
    dev.off()
    
    png("./Diablo Results/Correlation circle plot.png")
    plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
            pch = c(15,16), cex = c(2,2), col = c('darkorchid', 'lightgreen'))
    dev.off()
    
    #plotLoadings(sgccda.res, comp = 2, contrib = 'max', method = 'median')
    
    png("./Diablo Results/cim.png")
    cimDiablo(sgccda.res, size.legend = 0.5, legend.position = "topright", margins = c(10,15))
    dev.off()
    
    return(list(data.table1, data.table2))
    
  } else {
    print("Methyltion data introduced")
    
    #####################
    #STEP 1: PREPARATION#
    #####################
    
    # Check if the objects are data frames #
    
    if (class(Cancer_table) != "data.frame") {
      stop(print("Cancer_table is not a data.frame"))}
    
    if (class(mRNA) != "data.frame") {
      stop(print("mRNA is not a data.frame"))}
    
    if (class(CNV) != "data.frame") {
      stop(print("CNV is not a data.frame"))}
    
    if (class(Methylation) != "data.frame") {
      stop(print("Methylation is not a data.frame"))}
    
    
    # Set different gene names for each datatable #
    
    rownames(mRNA) <- paste(rownames(mRNA),"r",sep=".")
    rownames(CNV) <- paste(rownames(CNV),"c",sep=".") 
    rownames(Methylation) <- paste(rownames(Methylation),"m",sep=".") 
    
    
    ## Set samples names prepared for mixOmics with the same format and order ##
    
    # Remove totalmean column from Methylation table if it is present#
    
    if (colnames(Methylation)[1] == "totalmean"){
      Methylation <-Methylation[,-1]}
    
    # Change separations in samples names #
    
    colnames(mRNA) <- gsub("\\.", "-",colnames(mRNA))
    colnames(CNV) <- gsub("\\.","-",colnames(CNV))
    colnames(Methylation) <- gsub("\\.","-",colnames(Methylation))
    
    # Order to allow stage comparisons #
    
    mRNA <- mRNA[, order(match(colnames(mRNA), Cancer_table[,1]))]
    CNV <- CNV[, order(match(colnames(CNV), Cancer_table[,1]))]
    Methylation <- Methylation[, order(match(colnames(Methylation), Cancer_table[,1]))]
    
    if (all(sapply(list(colnames(mRNA),colnames(CNV),colnames(Methylation), Cancer_table[,1]), function(x) x == Cancer_table[,1])) == FALSE){
      stop(print("Samples are not the same"))}
    
    
    # Set inputs as numeric data frames #
    
    if (class(mRNA[1,1]) != c("numeric")){
      mRNA2 <- as.data.frame(sapply(mRNA, as.numeric))
      rownames(mRNA2) <- rownames(mRNA)   
      mRNA <- mRNA2
      print("Expression dataset changed to numeric data frame")}
    
    if (class(CNV[1,1]) != c("numeric")){
      CNV2 <- as.data.frame(sapply(CNV, as.numeric))
      rownames(CNV2) <- rownames(CNV)
      CNV <- CNV2
      print("CNV dataset changed to numeric data frame")}    
    
    if (class(Methylation[1,1]) != c("numeric")){
      Methylation2 <- as.data.frame(sapply(Methylation, as.numeric))
      rownames(Methylation2) <- rownames(Methylation)   
      Methylation <- Methylation2
      print("Methylation dataset changed to numeric data frame")}
    
    
    # Select the 10% of the genes with the highest standard deviation for each of the 3 data types #
    
    SD<-apply(mRNA,1,sd)
    top.10sd<-head(sort(SD,decreasing=TRUE),round(length(SD)/10))
    mRNA.f <- mRNA[names(top.10sd),]
    
    SD<-apply(CNV,1,sd)
    top.10sd<-head(sort(SD,decreasing=TRUE),round(length(SD)/10))
    CNV.f <- CNV[names(top.10sd),]
    
    SD<-apply(Methylation,1,sd)
    top.10sd<-head(sort(SD,decreasing=TRUE),round(length(SD)/10))
    Methylation.f <- Methylation[names(top.10sd),]
    
    
    #######################
    #STEP 2: APPLY  DIABLO#
    #######################
    
    
    
    # List of data sets #
    
    data = list(mRNA = t(mRNA.f),
                CNV = t(CNV.f),
                Methylation=t(Methylation.f))
    
    
    # Set Y to compare stage i vs stage iv #
    
    Y=as.factor(Cancer_table$tumor_stage)
    
    design = matrix(0.1, ncol = length(data), nrow = length(data), 
                    dimnames = list(names(data), names(data)))
    lapply(data,dim)
    Y
    diag(design) = 0 # To avoid analyzing an experiment against itself
    design
    
    
    sgccda.res = block.splsda(X = data, Y = Y, ncomp = 5, design = design) # Classifies a discrete outcome through the integration of 
    # multiple datasets selecting features from each data set, using a PLS-DA model
    
    
    # Evaluate performance of the PLS-DA model #
    
    set.seed(12345)
    
    perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)
    
    print("Saving results in Diablo Results folder in your working directory")
    dir.create(file.path(getwd(), "./Diablo Results"), showWarnings = FALSE)
    png("./Diablo Results/Diablo first model performance.png")
    plot(perf.diablo) # Components that exaplain most of the variability
    dev.off()
    
    ncomp = 2
    
    # Tune the keepX parameters obtained with the block.splsda function #
    
    s.keep <- sort(sample(1:length(rownames(Cancer_table)),5))
    test.keepX = list(mRNA = s.keep, CNV =s.keep, Methylation = s.keep)
    
    
    # Computes M-fold scores based on a user-input grid to determine the optimal parsity parameters #
    
    tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = 2,
                                  test.keepX = test.keepX, design = design,
                                  validation = 'Mfold', folds = 10, nrepeat = 1,
                                  cpus = 4, dist = "centroids.dist")
    
    list.keepX = tune.TCGA$choice.keepX
    
    
    # Apply method with tunned params #
    sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                              keepX = list.keepX, design = design)
    
    sgccda.res$design
    
    
    ################
    #STEP3: RESULTS#
    ################
    
    # Component 1 #
    
    varm1 <- head(selectVar(sgccda.res, comp = 1)$mRNA$value,100)
    varcnv1 <- head(selectVar(sgccda.res, comp = 1)$CNV$value,100)
    varmeth1 <- head(selectVar(sgccda.res, comp = 1)$meth$value,100)
    
    data.table1 <- rbind(varm1, varcnv1, varmeth1)
    data.table1 <- cbind(data.table1, rownames(data.table1))
    data.table1 <- data.table1[order(abs(data.table1$value.var), decreasing = TRUE),]
    colnames(data.table1) <- c("value.var1", "Genes")
    
    
    # Component 2 #
    
    varm2 <- head(selectVar(sgccda.res, comp = 2)$mRNA$value,100)
    varcnv2 <- head(selectVar(sgccda.res, comp = 2)$CNV$value,100)
    varmeth2 <- head(selectVar(sgccda.res, comp = 2)$meth$value,100)
    
    data.table2 <- rbind(varm2, varcnv2, varmeth2)
    data.table2 <- cbind(data.table2, rownames(data.table2))
    data.table2 <- data.table2[order(abs(data.table2$value.var), decreasing = TRUE),]
    colnames(data.table2) <- c("value.var2", "Genes")
    
    
    if (nrow(data.table1) > 100) {
      data.table1 <- data.table1[1:100,]}
    
    if (nrow(data.table2) > 100) {
      data.table2 <- data.table2[1:100,]}
    
    ###############
    #STEP4: PLOTS#
    ###############
    
    png("./Diablo Results/Diablo correlation first component.png")
    plotDiablo(sgccda.res, ncomp = 1)
    dev.off()
    
    png("./Diablo Results/Diablo correlation second component.png")
    plotDiablo(sgccda.res, ncomp = 2)
    dev.off()
    
    png("./Diablo Results/Arrow plot samples.png")
    plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
    dev.off()
    
    png("./Diablo Results/Correlation circle plot.png")
    plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
            pch = c(15,16, 17), cex = c(2,2,2), col = c('darkorchid', 'lightgreen',"red"))
    dev.off()
    
    #plotLoadings(sgccda.res, comp = 2, contrib = 'max', method = 'median')
    
    png("./Diablo Results/cim.png")
    cimDiablo(sgccda.res, size.legend = 0.5, legend.position = "topright", margins = c(10,15))
    dev.off()
    
    return(list(data.table1, data.table2))
    
  }
}

