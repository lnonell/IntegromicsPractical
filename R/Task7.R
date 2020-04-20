getwd()
# Task 7: Gerard Temprano Sagrera
# Purpose: Apply mixOmics to the data set comparing stage i vs stage iv.
# Input: Data-frames with mRNA, CN and methylation data and destination path. 
# Output: Inoformative graphs and the 100 most correlated variables.

task7<-function(Cancer_table, mRNA, CNV, Methylation, Destination_path){}

# Check if the objects are data frames:
if (class(Cancer_table) != "data.frame") {
  stop(print("Cancer_table is not a data.frame"))}

if (class(mRNA) != "data.frame") {
  stop(print("mRNA is not a data.frame"))}

if (class(CNV) != "data.frame") {
  stop(print("CNV is not a data.frame"))}

if (class(Methylation) != "data.frame") {
  stop(print("Methylation is not a data.frame"))}


# Check if the data frames contain the same samples:
#if(identical(rownames(mRNA), rownames(CNV)) == FALSE){
#  stop(print("The samples are not the same"))}

#if(identical(rownames(mRNA), rownames(Methylation)) == FALSE){
#  stop(print("The samples are not the same"))}

#if(identical(rownames(mRNA), rownames(Cancer_table)) == FALSE){
#  stop(print("The samples are not the same"))}

# Set different gene names for each datatable:
rownames(mRNA) <- paste(rownames(mRNA),"m",sep=".")
rownames(CNV) <- paste(rownames(CNV),"c",sep=".") 
rownames(Methylation) <- paste(rownames(Methylation),"me",sep=".") 

#Set samples names for mixOmics:
colnames(mRNA) <- sort(colnames(mRNA))
colnames(CNV) <- sort(colnames(CNV))
colnames(Methylation) <- sort(colnames(Methylation))

# Select the 10% of the genes with the highest standard deviation for each of the 3 data types:
SD<-apply(mRNA,1,sd)
top.10sd<-head(sort(SD,decreasing=TRUE),round(length(SD)/10))
mRNA.f <- mRNA[names(top.10sd),]

SD<-apply(CNV,1,sd)
top.10sd<-head(sort(SD,decreasing=TRUE),round(length(SD)/10))
CNV.f <- CNV[names(top.10sd),]

SD<-apply(Methylation,1,sd)
top.10sd<-head(sort(SD,decreasing=TRUE),round(length(SD)/10))
Methylation.f <- Methylation[names(top.10sd),]

# List of data sets:
data = list(mRNA = t(mRNA.f), 
            CNV = t(CNV.f),
            meth=t(Methylation.f))

lapply(data,dim)


# Set Y to compare stage i vs stage iv:
# Cancer_table
stage_iv <- grep("stage iv", LUSC.table$tumor_stage)
stage_i <- length(grep("stage i", LUSC.table$tumor_stage)) - length(stage_iv)

cond <- c(rep("stage i", stage_i), rep("stage iv", length(stage_iv)))
Y=as.factor(cond)

design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
lapply(data,dim)
Y
diag(design) = 0 # To avoid analyzing an experiment against itself
design


# Classify a discrete outcome through the integration of multiple datasets selecting features from each data set, using a PLS-DA model:
sgccda.res = block.splsda(X = data, Y = Y, ncomp = 5, design = design)

set.seed(123)
# Evaluate performance of the PLS-DA model:
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)

pdf(paste(Destination_path,"Diablo Performance 1st Model.pdf", sep = ""))  
pdf("Test.pdf")
plot(perf.diablo) # Components that exaplain most of the variability 
dev.off()

perf.diablo$choice.ncomp$WeightedVote
ncomp = 2

# Tune the keepX parameters obtained with the block.splsda function:
set.seed(1234) 
s.keep <- sort(sample(1:1770,7))
test.keepX = list(mRNA = s.keep, CNV =s.keep, Methylation = s.keep)

# Computes M-fold scores based on a user-input grid to determine the optimal parsity parameters 
tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp,
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              cpus = 4, dist = "centroids.dist")

list.keepX = tune.TCGA$choice.keepX


# Apply method with tunned params
sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design)

sgccda.res$design

###### Results
head(selectVar(sgccda.res)$mRNA$value,100)
head(selectVar(sgccda.res)$CNV$value,100)
head(selectVar(sgccda.res)$meth$value,100)

###### Plots
pdf(paste(Destination_path,"Diablo.pdf", sep = ""))  
pdf("Test.pdf")
plotDiablo(sgccda.res, ncomp = 1)

plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17), cex = c(2,2), col = c('darkorchid', 'lightgreen'))

plotLoadings(sgccda.res, comp = 2, contrib = 'max', method = 'median')

cimDiablo(sgccda.res, size.legend = 0.5, legend.position = "topright", margins = c(10,15))


dev.off()
