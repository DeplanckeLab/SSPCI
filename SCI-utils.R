
################## load prediction & aux files   #########################
load(paste(mpath,"/data/Prediction-Data1.RData",sep="")) #### load all data matrices (list) plus list of all gene intersection - allgenes & housekeeping genes - housekeep
load(paste(mpath,"/data/Prediction-Data2.RData",sep="")) #### load all data matrices (list) plus list of all gene intersection - allgenes & housekeeping genes - housekeep
load(paste(mpath,"/data/SSCI-TrainedModel.RData",sep=""))  #### trained model, as used in the publication
load(paste(mpath,"/data/HousekeepSet.RData",sep="")) #### set of housekeeping genes, as used in the publication
load(paste(mpath,"/data/Top3Models.RData",sep=""))  #### trained model, as used in the publication

nrq=20 ### initialise quantiles
mockcut=length(allgenes)/10 ### MAX NR of genes allowed missing in the matrix; adjust if required ####
choosecols=c("snow3","orange", "steelblue3"); mycex=0.5; choosecolsgenes=c("snow3","orange", "red") ### initialise colors ###

#### Filtering on the count matrix: remove low expr. genes & cells ######
#### take only subset of cells and genes from a dataset ####
standardSCFilterHousekeep=function(mmGenes,filgenes) {
        savegenes=mmGenes
        mmGenes=apply(mmGenes,2,function(x) as.numeric(as.vector(x))); rownames(mmGenes)=rownames(savegenes)
        in0=apply(mmGenes[which(rownames(mmGenes)%in%filgenes),],2,function(x) length(which(!x>0)))
	#### set the cutoff for column filtering: how many "housekeeping" genes are allowed to be absent; minimum present is 5 genes ####
	colcut=min(summary(in0)[5]+2*sd(in0),dim(mmGenes[which(rownames(mmGenes)%in%filgenes),])[1]-4)
        mmGenes=mmGenes[,which(apply(mmGenes[which(rownames(mmGenes)%in%filgenes),],2,function(x) length(which(!x>0)))<=colcut)]
        rowcut=unlist(mmGenes[which(rownames(mmGenes)%in%filgenes)[which(apply(mmGenes[which(rownames(mmGenes)%in%filgenes),],1,function(x) length(which(x>0)))/dim(mmGenes)[2]<=0.5)],])
        rowcut=rowcut[which(rowcut>0)]
        if (length(rowcut)%in%0) {
                rowcut=unlist(mmGenes[which(rownames(mmGenes)%in%filgenes)[which(apply(mmGenes[which(rownames(mmGenes)%in%filgenes),],1,function(x) length(which(x>0)))/dim(mmGenes)[2]<1)],])
                rowcut=summary(apply(rowcut,1,function(x) summary(x[which(x>0)])[2]))[2]
                } else {
                        rowcut=unlist(mmGenes[which(rownames(mmGenes)%in%filgenes)[which(apply(mmGenes[which(rownames(mmGenes)%in%filgenes),],1,function(x) length(which(x>0)))/dim(mmGenes)[2]<=0.5)],])
                        if(is.null(dim(rowcut))) {rowcut=summary(rowcut[which(rowcut>0)])[3]} else {
                                y=apply(rowcut,1,function(x) summary(x[which(x>0)])[5])
                                #rowcut=summary(rowcut)[5]
                                rowcut=summary(y)[3]}}
        mmGenes= mmGenes[apply(mmGenes, 1, function(x) length(x[x>rowcut])>=1),]
        return(mmGenes)
}

#### Wrap around prediction & plot generation  ####
##### Make predictions, save output prediction values to file & generate plot with prediction values #####
preppredict = function(predmatrix,testsetname,sepcats,mycols, rtsne_out=NULL,mytheta=0.00000001,myperp=40) {
	
	#### remove rows with NA or weird values #####
	predmatrix=predmatrix[which(!is.na(predmatrix[,1])),]
	predmatrix=predmatrix[which(!apply(predmatrix,1,function(x) sum(x))%in%c("NaN","-NaN","NA","Inf","-Inf")),]
	
	###### stop critical genes are not measured: l1.genes #####
	nomap.l1=length(which(!lasso.results$genes%in%rownames(predmatrix)))
	
	#### only continue if all critical genes present ####
	if (nomap.l1==0) {
		###### determine if there are IDs not present in the matrix  #####
		nomap=length(which(!allgenes%in%rownames(predmatrix)))
		if (nomap>0) {
			#### if the critical nr. of missing genes not reached, replace with mean values
			if ( nomap<=mockcut) {
				myval=apply(predmatrix,2,function(x) mean(x)) ##### motified this so it's not 0 [which(!x%in%0)]
				usematrix=rbind(predmatrix[allgenes[which(allgenes%in%rownames(predmatrix))],],
					t(sapply(1:nomap,function(x) myval)))
				rownames(usematrix)[which(!rownames(usematrix)%in%allgenes)]=allgenes[which(!allgenes%in%rownames(predmatrix))]
				usematrix=usematrix[allgenes,]	
			} else {usematrix=NULL}	
		} else {usematrix=predmatrix[allgenes,]}
	} else {usematrix=NULL}
	
	#### Generate predictions in case the matrix contains all required information 
	if (!is.null(usematrix)) {
		okmatrix.pred=apply(usematrix,2,function(x) { myq=c(0,quantile(x[which(!x%in%0)],  probs = c(1:nrq, NA)/nrq))
			z=sapply(x,function(y) max(which(y>=myq)));z[which(z%in%"-Inf")]=0; return(z)}) ### prepare matrix: transform to ventiles 

		#### Predict the new values, model is stored in lasso.results
		topred=as.matrix(t(okmatrix.pred)[,names(coef(lasso.results[[1]],s = "lambda.1se")[,1])[which(!names(coef(lasso.results[[1]],s = "lambda.1se")[,1])%in%"(Intercept)")]])
		mypred=predict(lasso.results[[1]],newx =topred,type="response",s = "lambda.1se") # Obtain "Stemness probability"
		mypred.cl=predict(lasso.results[[1]],newx = topred,type="class",s = "lambda.1se") # Obtain binary prediction	
		#### Write out predictions
		predout=cbind(rownames(mypred),mypred, mypred.cl); colnames(predout)=c("Cell","Stemness p","Binary SSC")
		write.table(predout,paste(mpath,"SSCI-predictions.",nrq,".",testsetname,".txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
		
		#### Get predictions per data subset: subcategories are determined by the "sepcats" variable
		remindex.dev=lapply(sepcats,function(x) colnames(okmatrix.pred[,grep(x,colnames(okmatrix.pred))]))
		
		#### Plot the results ####
		pdf(paste(mpath,"SSCI-predictions.",testsetname,".pdf",sep=""),width=8,height=8)
		par(mfrow=c(2,2))
		# box & barplots for the cell categories determines by sepcats
		remindex.dev.vals=lapply(remindex.dev,function(x) mypred[x,1])
		remindex.dev.cats=sapply(remindex.dev,function(x) c(length(which(mypred.cl[x,1]%in%"0")),
	 		length(which(mypred.cl[x,1]%in%"1"))))
		names(remindex.dev.vals)=colnames(remindex.dev.cats)=sepcats
		rownames(remindex.dev.cats)=c("0","1")
		boxplot(remindex.dev.vals,main=testsetname,col=c("orange"),las=2,cex.names=0.7,ylab="Stemness prob.")
		beanplot(remindex.dev.vals,main=testsetname,col=c("orange"),las=2,cex.names=0.7,what=c(1,1,1,0),
			overallline = "median",log="",bw="nrd0",ylab="Stemness prob.")
		barplot( apply(remindex.dev.cats,2,function(x) x/sum(x)),beside=F,main=testsetname,
			col=c("snow3",choosecols[3]),horiz=T,las=2,cex.names=0.7,xlab="Relative proportion")
		lapply(1:length(remindex.dev.vals),function(x) plot(density(remindex.dev.vals[[x]]),
			main=paste(testsetname,sepcats[x],sep=" "),xlim=c(0,1.5)))
		
		#### Determine, plot and write out top correlated genes with stemness probability
		temp=apply(predmatrix,1,function(x) cor(as.numeric(as.vector(mypred[colnames(predmatrix),1])),x,method="pearson"))
		temp.p=p.adjust(apply(predmatrix,1,function(x) cor.test(as.numeric(as.vector(mypred[colnames(predmatrix),1])),x,method="pearson")$p.val))
		temp=temp[which(temp.p<0.05)];temp.p=temp.p[which(temp.p<0.05)] #### top sign. correlated genes, at 5%
		if (length(temp[which(temp.p<0.05)])>0) {
			barplot(sort(temp,decreasing=T)[1:min(length(temp),10)],beside=T,las=2,cex.names=0.7,ylab="Pearson's cor",main=testsetname)
			write.table(sort(temp,decreasing=T),paste(mpath,"SSCI-predictions.",nrq,".",testsetname,"-correlations.txt",sep=""),quote=F,sep="\t")
		}
		
		#### 2D plots using either a provided or a newly calculated tSNE map; color according to cell categories, gene expression, or prediction results
		prefgenes=rownames(standardSCFilterHousekeep(predmatrix,c(housekeep.large,rownames(predmatrix)[grep("ERCC-",rownames(predmatrix))]))) ### susbet inf. genes
		if(length(prefgenes)%in%dim(predmatrix)[1]){prefgenes=rownames(predmatrix)[which(apply(predmatrix,1,function(x) length(which(x>=1)))>4)]}
		posmatrix=predmatrix[prefgenes,];mycols=mycols[colnames(posmatrix)]; mymax=max(posmatrix)*2/3
		rownames(posmatrix)=toupper(rownames(posmatrix)); set.seed=100
		mymatrix.withlab=posmatrix#[,which(allsums>1000)]
		if (!length(rtsne_out)>1) {
			rtsne_out <- Rtsne(as.matrix(t(scale(mymatrix.withlab))),theta=mytheta,perplexity=myperp) #### calculate RtSNE map, if not provided
		}
		rownames(rtsne_out$Y)=colnames(mymatrix.withlab); prefsne=rtsne_out; prefm=mymatrix.withlab
		# 2D plot; color according to cell categories
		plot(rtsne_out$Y,pch=19, main=paste(testsetname,sep=" "),col=mycols,xlab="Dimension 1",ylab="Dimension 2"); text(rtsne_out$Y, labels=colnames(mymatrix.withlab),col=mycols,cex=mycex)
		# color according to genes of interest
		interestg=unique(c(ginterest)); interestg=interestg[which(interestg%in%rownames(mymatrix.withlab))]
		for(i in 1:length(interestg)) {
			pcol=valuesToColorsAbs(unlist(mymatrix.withlab[interestg[i],]),choosecolsgenes,colnames(mymatrix.withlab),0,mymax); names(pcol)=colnames(mymatrix.withlab)
			plot(rtsne_out$Y,pch=19, main=paste(testsetname,interestg[i],sep=" "),col=pcol,xlab="Dimension 1",ylab="Dimension 2")
			text(rtsne_out$Y, labels=colnames(mymatrix.withlab),col=pcol,cex=mycex)
		}
		# 2D plot; color according to stemness probability
		plotcol=valuesToColorsAbs(unlist(mypred[,1][colnames(mymatrix.withlab)]),choosecols,colnames(mymatrix.withlab),0,1)
		plot(rtsne_out$Y,pch=19, main="Predicted",col=plotcol,xlab="Dimension 1",ylab="Dimension 2"); text(rtsne_out$Y, labels=colnames(mymatrix.withlab),
		xlab="Dimension 1",ylab="Dimension 2",col=plotcol,cex=mycex); legend("topleft",c("p=0","p=0.5","p=1"),col=choosecols,lwd=2)
		plotcol=c("snow3","steelblue3")[as.numeric(mypred.cl[colnames(mymatrix.withlab),1])+1]	
		# 2D plot; color according to binary stemness assignment
		plot(rtsne_out$Y,pch=19, main="Predicted-StandCut*",col=plotcol,xlab="Dimension 1",ylab="Dimension 2")
		text(rtsne_out$Y, labels=colnames(mymatrix.withlab),col=plotcol,cex=mycex); legend("topleft",c("p=0","p=1"),col=choosecols[c(1,3)],lwd=2)	
		dev.off()
		
		##### Return the RtSNE plot that was generated 
		return(rtsne_out)
	} else { return("Too many/too important genes missing in the input matrix; pls. make sure all classification genes are present & ~ all allgenes are there!")}
		
}

#### Wrap around prediction only  ####
#### Make predictions and return, no plotting #####
preppredict_simple = function(predmatrix,testsetname) {	
	
	#### remove rows with NA or weird values #####
	predmatrix=predmatrix[which(!is.na(predmatrix[,1])),]
	predmatrix=predmatrix[which(!apply(predmatrix,1,function(x) sum(x))%in%c("NaN","-NaN","NA","Inf","-Inf")),]
	
	###### stop critical genes are not measured: l1.genes #####
	nomap.l1=length(which(!lasso.results$genes%in%rownames(predmatrix)))
	
	#### only continue if all critical genes present ####
	if (nomap.l1==0) {
		###### determine if there are IDs not present in the matrix  #####
		nomap=length(which(!allgenes%in%rownames(predmatrix)))
		if (nomap>0) {
			#### if the critical nr. of missing genes not reached, replace with mean values
			if ( nomap<=mockcut) {
				myval=apply(predmatrix,2,function(x) mean(x)) ##### motified this so it's not 0 [which(!x%in%0)]
				usematrix=rbind(predmatrix[allgenes[which(allgenes%in%rownames(predmatrix))],],
					t(sapply(1:nomap,function(x) myval)))
				rownames(usematrix)[which(!rownames(usematrix)%in%allgenes)]=allgenes[which(!allgenes%in%rownames(predmatrix))]
				usematrix=usematrix[allgenes,]	
			} else {usematrix=NULL}	
		} else {usematrix=predmatrix[allgenes,]}
	} else {usematrix=NULL}
	
	#### Generate predictions in case the matrix contains all required information 
	if (!is.null(usematrix)) {
		okmatrix.pred=apply(usematrix,2,function(x) { myq=c(0,quantile(x[which(!x%in%0)],  probs = c(1:nrq, NA)/nrq))
			z=sapply(x,function(y) max(which(y>=myq)));z[which(z%in%"-Inf")]=0; return(z)}) ### prepare matrix: transform to ventiles 

		#### Predict the new values, model is stored in lasso.results
		topred=as.matrix(t(okmatrix.pred)[,names(coef(lasso.results[[1]],s = "lambda.1se")[,1])[which(!names(coef(lasso.results[[1]],s = "lambda.1se")[,1])%in%"(Intercept)")]])
		mypred=predict(lasso.results[[1]],newx =topred,type="response",s = "lambda.1se") # Obtain "Stemness probability"
		mypred.cl=predict(lasso.results[[1]],newx = topred,type="class",s = "lambda.1se") # Obtain binary prediction	
		#### Write out predictions
		predout=cbind(rownames(mypred),mypred, mypred.cl); colnames(predout)=c("Cell","Stemness p","Binary SSC")
		write.table(predout,paste(mpath,"SSCI-predictions.",nrq,".",testsetname,".txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
		return(predout)
	} else { return("Too many/too important genes missing in the input matrix; pls. make sure all classification genes are present & ~ all allgenes are there!")}
	
}

#### Wrap around prediction only - for top three best performing models (lasso, elastic-net, randomforest & final 23 gene model) ####
#### Make predictions and return, no plotting #####
#### top3models needs to be loaded - contains the models ######
preppredict_models=function(predmatrix,testsetname) {
	#### remove rows with NA or weird values #####
	predmatrix=predmatrix[which(!is.na(predmatrix[,1])),]
	predmatrix=predmatrix[which(!apply(predmatrix,1,function(x) sum(x))%in%c("NaN","-NaN","NA","Inf","-Inf")),]

	###### stop critical genes are not measured: l1.genes #####
	nomap.l1=length(which(!lasso.results$genes%in%rownames(predmatrix)))

	#### only continue if all critical genes present ####
	if (nomap.l1==0) {
		###### determine if there are IDs not present in the matrix  #####
		nomap=length(which(!allgenes%in%rownames(predmatrix)))
		if (nomap>0) {
			#### if the critical nr. of missing genes not reached, replace with mean values
			if ( nomap<=mockcut) {
				myval=apply(predmatrix,2,function(x) mean(x)) ##### motified this so it's not 0 [which(!x%in%0)]
				usematrix=rbind(predmatrix[allgenes[which(allgenes%in%rownames(predmatrix))],],
				t(sapply(1:nomap,function(x) myval)))
				rownames(usematrix)[which(!rownames(usematrix)%in%allgenes)]=allgenes[which(!allgenes%in%rownames(predmatrix))]
				usematrix=usematrix[allgenes,]
				} else {usematrix=NULL}
				} else {usematrix=predmatrix[allgenes,]}
				} else {usematrix=NULL}

#### Generate predictions in case the matrix contains all required information 
if (!is.null(usematrix)) {
	okmatrix.pred=apply(usematrix,2,function(x) { myq=c(0,quantile(x[which(!x%in%0)],  probs = c(1:nrq, NA)/nrq))
		z=sapply(x,function(y) max(which(y>=myq)));z[which(z%in%"-Inf")]=0; return(z)}) ### prepare matrix: transform to ventiles 

		topred=t(okmatrix.pred)
		lasso.pred.p=predict(top3models[[1]],newx=topred,type="response",s = "lambda.1se")
		lasso.pred.bin=predict(top3models[[1]],newx=topred,type="class")
		elastic.pred.p=predict(top3models[[2]],newx=topred,type="response",s = "lambda.1se")
		elastic.pred.bin=predict(top3models[[2]],newx=topred,type="class")
		random.pred.p=predict(top3models[[3]],newdata=topred,type="prob")[,2]
		random.pred.bin=predict(top3models[[3]],newdata=topred,type="class")
		myprobs=cbind(lasso.pred.p,lasso.pred.bin,elastic.pred.p,elastic.pred.bin,random.pred.p,as.numeric(random.pred.bin)-1)

		#### Predict the new values, model is stored in lasso.results
		topred=as.matrix(t(okmatrix.pred)[,names(coef(lasso.results[[1]],s = "lambda.1se")[,1])[which(!names(coef(lasso.results[[1]],s = "lambda.1se")[,1])%in%"(Intercept)")]])
		mypred=predict(lasso.results[[1]],newx =topred,type="response",s = "lambda.1se") # Obtain "Stemness probability"
		mypred.cl=predict(lasso.results[[1]],newx = topred,type="class",s = "lambda.1se") # Obtain binary prediction

		#### Write out predictions
		predout=cbind(rownames(mypred),myprobs,mypred,mypred.cl); colnames(predout)=c("Cell","Lasso.p","Lasso.bin",
		"Elastic.p","Elastic.bin","Random.p","Random.bin","Lasso.23.p","Lasso.23.bin")
		write.table(predout,paste(mpath,"SSPCI-models-predictions.",nrq,".",testsetname,".txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
		return(predout)
		} else { return("Too many/too important genes missing in the input matrix; pls. make sure all classification genes are present & ~ all allgenes are there!")}
	}
	
#### Expression values to colors for tSNE plots ####
valuesToColorsAbs=function (myvector,choosecols,mynames,mymin,mymax) {
	myvector=c(mymin,mymax,myvector); names(myvector)[1:2]=c("A1","A2"); mynames=c("A1","A2",mynames)
	myvector[which(is.na(as.numeric(as.vector(myvector))))]=min(myvector[which(!is.na(as.numeric(as.vector(myvector))))])
	myvector=round(myvector,digits=2)
	colormap=seq(min(myvector),max(myvector),by=0.01)
	colormap=colorRampPalette(choosecols)(length(colormap))
	names(colormap)=round(seq(min(myvector),max(myvector),by=0.01),digits=2)
	rpfmycol=as.vector(colormap[match(myvector,names(colormap))])
	names(rpfmycol)=mynames
	return(rpfmycol[which(!names(rpfmycol)%in%c("A1","A2"))])
}

##### Read the count matrix stored in a .csv file, htseq input #####
##### Output a list with counts/gene, colors & cell cycle markers ############# 
readFile = function(filename,inputtype="simple",sublist=NULL) {
	filecols=brewer.pal(9,"Pastel1") #### change this for different colors 
	leng=read.csv(paste(mpath,filename,sep=""))
	### simple input matrix is default, but read like this if htseq counts format ###
	if (inputtype%in%"htseq") {
		leng=leng[which(!leng[,1]%in%c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")),]
		leng.spike=leng[grep("ERCC-",leng[,1]),]; rownames(leng.spike)=leng.spike[,1] #### if there are spike-ins
		leng.gene=leng[which(!leng[,1]%in%leng.spike[,1]),]
		leng.gene=leng.gene[which(!duplicated(toupper(leng.gene[,1]))),]
		rownames(leng.gene)=toupper(leng.gene[,1])
		leng=rbind(leng.gene,leng.spike); leng=leng[,-1]
		myrownames=rownames(leng); mycolnames=colnames(leng)
		leng=apply(leng,2,function(x) as.numeric(unlist(x)))
		colnames(leng)=mycolnames;rownames(leng)=myrownames      
		} 	
	##### list of potentially interesting genes - here cell cycle - change if others are interesting #####
	mymarker=toupper(c("Mki67","Pdpn","Aurkb","Ccnd1","Ccnd2","Cdk2","Cdk4","Cep110","Cep250","E2f2","E2f3","E2f7","E2f8",
		"Mapk1","Mapk12","Mapk13","Mapk3","Mapk4","Mapk7","Mcm2","Mcm3","Mcm4","Mcm6","Src"))
	##### set colors ######
	mycol=rep("snow3",length(colnames(leng)))
	if (length(sublist)>0) {
		lengcol=sapply(sublist,function(x) grep(x,colnames(leng)))
		for (x in 1:length(lengcol)) {mycol[lengcol[[x]]]=filecols[x]}
	} 
	names(mycol)=colnames(leng);lenglist=list(leng,mycol,mymarker)
	names(lenglist)=c("leng","Colors","Marker"); return(lenglist)
}	

##### Read the pancreas Baaron data #####
##### Output a list with counts/gene & categories ############# 
read.baaron=function(myfile,myname) {
	pcn1=t(read.table(paste(mpath,myfile,sep=""),head=T,sep=","))
	colnames(pcn1)=paste(myname,mysplit(mysplit(pcn1[1,],"_",4,2),"[.]",2,1),"C",mysplit(pcn1[1,],"_",4,4),sep="")
	mycats=pcn1[3,]; pcn1=pcn1[-c(1:3,dim(pcn1)[1]),]; rownames(pcn1)=toupper(rownames(pcn1))
	savedrow=rownames(pcn1);savedcol=colnames(pcn1);
	pcn1=apply(pcn1,2,function(x) as.numeric(as.vector(x)))
	rownames(pcn1)=savedrow
	mylist=list(pcn1,mycats); names(mylist)=c("Matrix","Cats")
	return(mylist)
}

##### Auxiliary function, used for spplitting a vector of strings ######
#### e.g. me-2 mB-3 tip-4 will be split mysplit(vector,"-",2,2) returns 2 3 4
mysplit <- function (vector,splitsign,numbers,position) {
	return(unlist(strsplit(as.character(vector),splitsign))[seq(position,length(vector)*numbers,by=numbers)])
}
