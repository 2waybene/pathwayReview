#source("http://bioconductor.org/biocLite.R")
#biocLite("CePa")
library("ROntoTools")
#library(CePa)
#data(PID.db)
#library(org.Hs.eg.db)
path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
#setwd("/wsu/home/gd/gd03/gd0393/Pathway/dataset/");
#path="/wsu/home/gd/gd03/gd0393/Pathway/dataset/"


#dataSets=c("GSE1297","GSE28146","GSE5281_EC","GSE5281_HIP","GSE5281_MTG","GSE5281_PC", "GSE5281_SFG", "GSE5281_VCX")
#load("diseaseList.RData")
#dataSets=c("GSE1297","GSE28146","GSE5281_EC","GSE5281_HIP","GSE5281_MTG","GSE5281_PC", "GSE5281_SFG", "GSE5281_VCX")
dataSets=c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
           "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")
#dataSets <- dataSets[Alzheimers_disease]
#dataSets=c("GSE781")
iterno=1000
PVAL <- 0.05
maxDE=400
load("KEGG65.150Pathways.RData")
#load("/wsu/home/gd/gd03/gd0393/Pathway/KEGG65.150Pathways.RData")
kpg <- setEdgeWeights(kpg)
kpg <- setNodeWeights(kpg, defaultWeight = 1)

getnodes=function(x){
	return (x@nodes);}

keggnodes=lapply(kpg,FUN=function(x){return (x@nodes)})
allGenes=unique(unlist(keggnodes))
#mymapping=as.data.frame(org.Hs.egSYMBOL)
#mygeneSym=as.character(mymapping$symbol)
#names(mygeneSym)=mymapping$gene_id
#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/mygeneSym.RData")
#load("/wsu/home/gd/gd03/gd0393/Pathway/mygeneSym.RData")


myData <- list()
for(i in 1:(length(dataSets))) {
  
	dataset=dataSets[i]
	print(dataset)
	load(paste(path,"dataset/", dataset,"/",dataset,".RData",sep=""))
	#load(paste(path, dataset,"/",dataset,".RData",sep=""))
	group <- get(paste("group_",dataset,sep=""))
	annotation=get(paste("annotation_",dataset,sep=""))
  	#design=get(paste("design_",dataset,sep=""))
	design=c("Not Paired")

	ano=paste(substring(annotation, 1, nchar(annotation)-3),"ENTREZID",sep="")
	require(annotation, character.only = TRUE)
	geneAno=as.data.frame(get(ano))
	geneAnno=geneAno[,2]
	names(geneAnno)=geneAno[,1]
	allProbes=unique(geneAno$probe_id)

	data <- 2^get(paste("data_",dataset,sep=""))
	data <- data[rownames(data)%in%allProbes,]
	rownames(data) <- as.character(geneAnno[rownames(data)])
	mydat <- aggregate(data,by=list(rownames(data)),FUN=median)
	mydat <- mydat[paste("hsa:",mydat$Group.1,sep="")%in%allGenes,]


	controlDat = mydat[,rownames(group[group$Group=="c",])]
	rownames(controlDat) = as.character(mydat$Group.1)

	diseaseDat = mydat[,rownames(group[group$Group=="d",])]
	rownames(diseaseDat) = as.character(mydat$Group.1)


	controlMean <- log(apply(controlDat,MARGIN=1,mean),base=2)
	diseaseMean <- log(apply(diseaseDat,MARGIN=1,mean),base=2)

	foldChange <- diseaseMean-controlMean
	names(foldChange) <- paste("hsa:",names(foldChange),sep="")



	#if paired, needs some attention
	pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
	names(pvalues) <- names(foldChange)

	#DEGenes <- foldChange[(abs(foldChange)>FC) & (pvalues<PVAL)]
	DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
  #maxDE=min(floor(length(allGenes)/10),length(DEGenes))
  DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
  DEGenes <- foldChange[names(DEGenes)]
  
  temp=data.frame(row.names=names(kpn), Name=kpn,pvalue=rep(NA,length(kpn)),rank=rep(NA,length(kpn)))
  for (pw in 1:length(kpn)) {
    #print(kpn[pw])
    #pwName <- names(kpn)[pw]
    #DEhit <- DEGenes[intersect(get(paste0("keggnodes$`",pwName,"`")),names(DEGenes))]
    DEhit <- DEGenes[intersect(keggnodes[[pw]],names(DEGenes))]
    DEmiss <- DEGenes[!DEGenes%in%DEhit]
    
    
    if (length(DEhit) == 0 || length(DEmiss) == 0) { pval = 1;}
    else {
      ks <- ks.test(DEhit,DEmiss)
      pval = ks$p.value
    }
  
    temp[names(kpn)[pw],"pvalue"] <- pval
  }
  temp$rank <- rank(temp[,2],ties.method="min")
  myData[[i]] <- temp
}

names(myData) <- dataSets
save(myData,file=paste(path,"alldatasetsKS.RData",sep=""))

################ number of significant pathway ###########

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsKS.RData")
signPW = vector(mode = "numeric", length = 0)
for (i in 1:length(myData)) {
  signPW = c(signPW, length(which(myData[[i]]$pvalue<0.05))/150)
}

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/KS40Datasets.RData")
for (i in 1:length(myData)) {
  signPW = c(signPW, length(which(myData[[i]]$pvalue<0.05))/150)
}

mean(signPW)
