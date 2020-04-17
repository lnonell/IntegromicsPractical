#Task 5 Óscar Baeza
#Purpose: Perform correlations (Spearman) (CN/mRNA and mRNA/methylation) and retrieve a table with good correlated genes and regression plots.
#Input: 2 data.frames with format being samples in columns and genes (HUGO Symbol) in rows (mRNA data),
# samples in columns and genes (HUGO Symbol) in rows containing the segment mean of each gene and sample (CN data),
# samples in columns and genes (HUGO Symbol) in rows containing the beta values mean of each gene and sample (methilation data).
#Output: data.frame with good correlations between data.frames (in absolute value >0.67)
task5<-function(dataframe1, dataframe2, method){
  return(correlations)
}

dataframe1 <- data.frame(mtcars)
dataframe1 <- c("1","2")

if (class(dataframe1) == "data.frame") {
  print("It is a dataframe")
  } else {
    stop(paste(deparse(substitute(dataframe1)))," is not a data.frame")
}


method <- "CN/mRNA"
if(method == "CN/mRNA") {
  print("Hello")
}

KEKW