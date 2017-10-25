################## Load the trained model & provides examples of how to generate predictions on 4 datasets  #########################

###### load required libraries, paths & functions ######
###### If not installed, install:
#install.packages("glmnet"); install.packages("randomForest")
library(glmnet); library(randomForest)
#mpath=""; ##### set the path, otherwise default
source(paste(mpath,"/SCI-utils.R",sep="")) ##### read the necessary functions
mpath=paste(mpath,"/data/",sep="") ##### change read/write folder to mpath/data

##### EXAMPLE.1: Schlitzer et al. 2015 (dendritic cell differentiation)
testsetname="SchlitzAll"
mymatrix=standardSCFilterHousekeep(schlitzlist[[1]],c(housekeep.large)) #### filter cells  with low counts overall, replace or remove depending on matrix  
mymatrix=log1p(schlitzlist[[1]][unique(c(allgenes[which(allgenes%in%rownames(schlitzlist[[1]]))],rownames(mymatrix))),colnames(mymatrix)])
schlitzpredict=preppredict_simple(mymatrix,paste(testsetname,"-SSCI",sep="")) #### generate predictions
schlitzpredict.allmodels=preppredict_models(mymatrix,paste(testsetname,"-SSCI",sep="")) #### generate using lasso & elastic-net on all genes; randomForest on all genes & SSPCI 

##### EXAMPLE.2: . New matrix: read from file; Leng et al. 2015, H1.ESC cells & cells in different cell cycle stage 
testsetname="LengAll"
lengList=readFile("GSE64016_H1andFUCCI_normalized_EC.csv","htseq")
mymatrix=standardSCFilterHousekeep(lengList[[1]],c(housekeep.large)) #### filter cells with low counts overall, replace or remove depending on matrix  
mymatrix=log1p(lengList[[1]][unique(c(allgenes[which(allgenes%in%rownames(lengList[[1]]))],rownames(mymatrix))),colnames(mymatrix)])
lengpredict=preppredict_simple(mymatrix,paste(testsetname,"-SSCI",sep="")) #### generate predictions
lengpredict.allmodels=preppredict_models(mymatrix,paste(testsetname,"-SSCI",sep="")) #### generate using lasso & elastic-net on all genes; randomForest on all genes & SSPCI 
