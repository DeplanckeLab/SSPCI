################## Load the trained model & provides examples of how to generate predictions on 4 datasets  #########################

###### load required libraries ######
###### If not installed, install:
#install.packages("glmnet"); install.packages("Rtsne"); install.packages("randomForest"); install.packages("beanplot"); install.packages("RColorBrewer") 
library(glmnet); library(Rtsne);  library(randomForest); library(beanplot); library(RColorBrewer)

#mpath=""; # set the path, otherwise default #####
source(paste(mpath,"SCI-utils.R",sep=""))
##### change read/write folder
mpath=paste(mpath,"data/",sep="")

################## load prediction & tSNE maps used in the publication   #########################
#load(paste(mpath,"Prediction-matrices.RData",sep="")) ###### matrices the model was trained on, as used in the publication; only load if required
#load(paste(mpath,"tSNE-PlotsMatrices.RData",sep="")) #### tSNE plots & matrices used in the publication; only load if required

################## Initialise some variables ##########
choosecols=c("snow3","orange", "steelblue3"); mycex=0.5; choosecolsgenes=c("snow3","orange", "red") ### initialise colors ###
ginterest=sort(toupper(c("Lgr5","Mki67",lasso.results$genes))) ### Genes whose expression will be plotted ###

##### EXAMPLE.1: Schlitzer et al. 2015 (dendritic cell differentiation)
testsetname="SchlitzAll"
mymatrix=standardSCFilterHousekeep(schlitzlist[[1]],c(housekeep.large)) #### filter cells  with low counts overall, replace or remove depending on matrix  
mymatrix=log1p(schlitzlist[[1]][unique(c(allgenes[which(allgenes%in%rownames(schlitzlist[[1]]))],rownames(mymatrix))),colnames(mymatrix)])
sepcats=c("MDC","CDP","PreDC"); mycols=schlitzlist[["Colors"]] #### determine the three categories to separately analyse
schlitzpredict=preppredict_simple(mymatrix,paste(testsetname,"-SSCI",sep="")) #### generate predictions only
schlitzpredict.allmodels=preppredict_models(mymatrix,paste(testsetname,"-SSCI",sep="")) #### generate using lasso & elastic-net on all genes; randomForest on all genes & SSPCI 
schltestsne=preppredict(mymatrix,paste(testsetname,"-SSCI",sep=""),sepcats,mycols) #### generate predictions & plots as in manuscript
save(schltestsne,file=paste(mpath,"SSCI-",testsetname,".Predictions.RData",sep="")) ### save tSNE to keep the same map

##### EXAMPLE.2:  Kowalczyk et al. 2015; this dataset is very large => use only subset of genes
testsetname="KowaSub"
mymatrix=standardSCFilterHousekeep(kowalist[[1]],c(housekeep.large)) #### filter cells with low counts overall, replace or remove depending on matrix  
mymatrix=kowalist[[1]][unique(c(allgenes[which(allgenes%in%rownames(kowalist[[1]]))],rownames(mymatrix))),colnames(mymatrix)]
##### Subsample data ####
if (dim(mymatrix)[2]>500) {
	okmatrixCut=summary(apply(mymatrix,2,function(x) length(which(x>1)))); okmatrix=mymatrix
	okmatrix=okmatrix[,which(apply(okmatrix,2,function(x) length(which(x>1)))>=okmatrixCut[2])]
	okmatrix=okmatrix[,which(apply(okmatrix,2,function(x) length(which(x>1)))<=okmatrixCut[5])]
	ind1=colnames(okmatrix); ind1=ind1[seq(1,length(ind1),by=round(length(ind1)/300))]
	mymatrix=okmatrix[,c(ind1)]
}
sepcats=c("yST.HSC","oST.HSC","yLT.HSC","oLT.HSC","y.MPP","o.MPP"); mycols=kowalist[["Colors"]][colnames(mymatrix)] #### set colors and categories
kowapredict=preppredict_simple(mymatrix,paste(testsetname,"-SSCI",sep="")) #### generate predictions only
kowapredict.allmodels=preppredict_models(mymatrix,paste(testsetname,"-SSCI",sep="")) #### generate using lasso & elastic-net on all genes; randomForest on all genes & SSPCI 
kowatestsne=preppredict(mymatrix,paste(testsetname,"-SSCI",sep=""),sepcats,mycols,NULL,0.00000001,80) #### generate predictions & plots as in manuscript
save(kowatestsne,file=paste(mpath,"SSCI-",testsetname,".Predictions.RData",sep="")) ### save tSNE to keep the same map

##### EXAMPLE.3: . New matrix: read from file; Leng et al. 2015, H1.ESC cells & cells in different cell cycle stage 
sepcats=c("H1_","G2_","S_","G1_"); testsetname="LengAll"
lengList=readFile("GSE64016_H1andFUCCI_normalized_EC.csv","htseq",sepcats)
mymatrix=standardSCFilterHousekeep(lengList[[1]],c(housekeep.large)) #### filter cells with low counts overall, replace or remove depending on matrix  
mymatrix=log1p(lengList[[1]][unique(c(allgenes[which(allgenes%in%rownames(lengList[[1]]))],rownames(mymatrix))),colnames(mymatrix)])
mycols=lengList[["Colors"]] #### use the colors 
lengpredict=preppredict_simple(mymatrix,paste(testsetname,"-SSCI",sep="")) #### generate predictions only
lengpredict.allmodels=preppredict_models(mymatrix,paste(testsetname,"-SSCI",sep="")) #### generate using lasso & elastic-net on all genes; randomForest on all genes & SSPCI 
lengtestsne=preppredict(mymatrix,paste(testsetname,"-SSCI",sep=""),sepcats,mycols)  #### generate predictions & plots as in manuscript
save(lengtestsne,file=paste(mpath,"SSCI-",testsetname,".Predictions.RData",sep=""))

