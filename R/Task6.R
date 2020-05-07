#Task 6: Blanca Rodriguez Fernandez 
#Purpose: Apply MFA to the mRNA, CN and methylation data comparing stage iv vs stage i.
#input: files of task1, task2, task3, task4 + path for output files
#outputs: dataframe 100 most correlated variables with PC1 and PC2 + MFA plots in png format

task6 <- function(df_samples, df.rna, df.cn, df.met, pth = getwd(),...){
 
  ###########################
  ## Load needed libraries ##
  ###########################
  
  suppressPackageStartupMessages(library(FactoMineR))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(factoextra))
  suppressPackageStartupMessages(library(stringr))
  
  ##################################
  ## Define filter by SD function ##
  ##################################
  filterSD <- function(data, percentage = 0.1){
    SD <- apply(data,1,sd)
    top10sd <- head(sort(SD,decreasing = TRUE), round(nrow(data)*percentage))
    data.f <- data[names(top10sd),]
    return(data.f)
  }
  
  ##################################
  ## Transform CN to class matrix ##
  ##################################
  n.cn <- apply(df.cn,2, as.numeric) #class numeric is needed
  rownames(n.cn) <- rownames(df.cn)
    
  ## Now, we can perfom MFA w/o methylation data. If methylation data  
  ## is missing, MFA will be applied to CN and expression data. 
  
  if(missing(df.met)){
    
    ##########################
    ## MFA without methylation
    ##########################
    
    warning(print("Methylation data is missing"))
    
    if (pth == getwd()){
      warning(print("Default output file is your current working directory"))
    }
    
    ## Check arguments ##
    #####################
    stopifnot(is.data.frame(df_samples))
    stopifnot(is.data.frame(df.rna))
    stopifnot(is.data.frame(df.cn))
    stopifnot(is.numeric(n.cn[,1]))
    
    
    ## Filter 10% genes by standard deviation ##
    ############################################
    rna.f <- filterSD(df.rna)
    cn.f <- filterSD(n.cn)
    
    ## Set colnames order equal to task1 ##
    #######################################
    rna.f <- rna.f[, order(match(colnames(rna.f), as.character(df_samples[,1])))]
    cn.f <- cn.f[, order(match(colnames(cn.f), as.character(df_samples[,1])))]
    
    if(identical(colnames(rna.f), colnames(cn.f)) != TRUE){
      stop(print("Samples are not the same"))
    }
    
    ## Data preparation for MFA ##
    ##############################
    ## Remove NAs 
    rna4MFA <- rna.f[!is.na(rna.f[,1]),]
    cn4MFA <- cn.f[!is.na(cn.f[,1]),]
    
    ## Label genes: c = copy number; r = expression data
    rownames(rna4MFA) <- paste(rownames(rna4MFA),"r",sep=".")
    rownames(cn4MFA) <- paste(rownames(cn4MFA),"c",sep=".") 
    
    ## Define conditions 
    cond <- as.factor(df_samples$tumor_stage)
    
    ## Number of genes in each dataset 
    rna.l<-nrow(rna4MFA)
    cn.l<-nrow(cn4MFA)
    
    ## New data frame with individuals in rows and variables in columns
    data4Facto<-data.frame(cond,t(rna4MFA),t(cn4MFA)) 
    cat("\nGetting participant identifier\n")
    rownames(data4Facto) <- paste(str_sub(df_samples$barcode,-4), cond, sep="_")
    
    ## Apply MFA ##
    ###############
    res.cond <- MFA(data4Facto, group=c(1,rna.l,cn.l), type=c("n","c","c"), 
                    ncp=5, name.group=c("cond","RNA","CN"),num.group.sup=c(1), graph = FALSE) 
    
    ## Obtain informative plots ##
    ##############################
    # Group of variables
    fviz_mfa_var(res.cond, "group")
    ggsave(file=paste(pth,"VariableGroupsMFA.png",sep ="/"))
    
    # Partial individuals
    fviz_mfa_ind(res.cond,habillage = cond, palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                 addEllipses = TRUE, ellipse.type = "confidence", 
                 repel = TRUE 
    ) 
    ggsave(file=paste(pth,"IndividualsMFA.png",sep ="/"))
    
    # Partial axes
    fviz_mfa_axes(res.cond, palette = c("#00AFBB", "#E7B800", "#FC4E07"))
    ggsave(file=paste(pth,"PartialAxesMFA.png",sep ="/"))
    
    cat("\nCheck MFA plots in your working directory or output path\n")
    
    ## Return 100 most correlated variables ##
    ##########################################
    ## PC1 (dimension 1 of global PCA)
    PC1 <- names(head(sort(res.cond$global.pca$var$cor[,1],decreasing = TRUE),100))
    ## PC2 (dimension 2 of global PCA)
    PC2 <- names(head(sort(res.cond$global.pca$var$cor[,2],decreasing = TRUE), 100))
    cat("\nThese are the 100 most correlated genes to PC1 and PC2\n")
    return(data.frame(PC1,PC2)) 
    
  } else{
    
    ##########################
    ## MFA with methylation ##
    ##########################

    ## Check arguments ##
    #####################
    stopifnot(is.data.frame(df_samples))
    stopifnot(is.data.frame(df.rna))
    stopifnot(is.data.frame(df.cn))
    stopifnot(is.numeric(n.cn[,1]))
    stopifnot(is.numeric(df.met[,1]))
    
    ## Drop total mean column from methylation data ##
    ##################################################
    df.met <- df.met[,-1]
    
    ## Filter 10% genes by standard deviation ##
    ############################################
    rna.f <- filterSD(df.rna)
    cn.f <- filterSD(n.cn)
    met.f <- filterSD(df.met)
    
    ## Set colnames order equal to task1 ##
    #######################################
    rna.f <- rna.f[, order(match(colnames(rna.f), as.character(df_samples[,1])))]
    cn.f <- cn.f[, order(match(colnames(cn.f), as.character(df_samples[,1])))]
    met.f <- met.f[, order(match(colnames(met.f), as.character(df_samples[,1])))]
    
    if(identical(colnames(met.f),identical(colnames(rna.f), colnames(cn.f))) != TRUE){
      stop(print("Samples are not the same"))
    }
    
    ## Data preparation for MFA ##
    ##############################
    ## Remove NAs 
    rna4MFA <- rna.f[!is.na(rna.f[,1]),]
    cn4MFA <- cn.f[!is.na(cn.f[,1]),]
    met4MFA <- met.f[!is.na(met.f[,1]),]
    
    ## Label genes: c = copy number; r = expression data; m = methylation
    rownames(rna4MFA) <- paste(rownames(rna4MFA),"r",sep=".")
    rownames(cn4MFA) <- paste(rownames(cn4MFA),"c",sep=".") 
    rownames(met4MFA) <- paste(rownames(met4MFA),"m",sep=".")
    
    ## Define conditions 
    cond <- as.factor(df_samples$tumor_stage)
    
    ## Number of genes in each dataset 
    rna.l<-nrow(rna4MFA)
    cn.l<-nrow(cn4MFA)
    met.l<-nrow(met4MFA)
    
    # New data frame with individuals in rows and variables in columns
    data4Facto<-data.frame(as.factor(cond),t(rna4MFA),t(cn4MFA),t(met4MFA)) 
    cat("\nGetting participant identifier\n")
    rownames(data4Facto) <- paste(str_sub(df_samples$barcode,-4), cond, sep="_")
    
    ## Apply MFA ##
    ###############
    res.cond <- MFA(data4Facto, group=c(1,rna.l,cn.l,met.l), type=c("n","c","c","c"), 
                    ncp=5, name.group=c("cond","RNA","CN","MET"),num.group.sup=c(1), graph = FALSE) 
    
    ## Obtain informative plots ##
    ##############################
    # Group of variables
    fviz_mfa_var(res.cond, "group")
    ggsave(file=paste(pth,"VariableGroupsMFA.png",sep ="/"))
    
    # Partial individuals
    fviz_mfa_ind(res.cond,habillage = cond, palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                 addEllipses = TRUE, ellipse.type = "confidence", 
                 repel = TRUE 
    ) 
    ggsave(file=paste(pth,"IndividualsMFA.png",sep ="/"))
    
    # Partial axes
    fviz_mfa_axes(res.cond, palette = c("#00AFBB", "#E7B800", "#FC4E07"))
    ggsave(file=paste(pth,"PartialAxesMFA.png",sep ="/"))
    cat("\nCheck MFA plots in your working directory or output path\n")
    ## Return 100 most correlated variables with 
    ## PC1 (dimension 1 of global PCA)
    PC1 <- names(head(sort(res.cond$global.pca$var$cor[,1],decreasing = TRUE),100))
    ## PC2 (dimension 2 of global PCA)
    PC2 <- names(head(sort(res.cond$global.pca$var$cor[,2],decreasing = TRUE), 100))
    cat("\nThese are the 100 most correlated genes to PC1 and PC2\n")
    return(data.frame(PC1,PC2))
    
  }
    
}

