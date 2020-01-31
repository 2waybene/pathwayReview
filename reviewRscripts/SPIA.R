setwd("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/")
path="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/"
# setwd("/wsu/home/gd/gd03/gd0393/Pathway/dataset/");
# path="/wsu/home/gd/gd03/gd0393/Pathway/dataset/"

iterno=2000
PVAL <- 0.5
maxDE=400

source("https://bioconductor.org/biocLite.R")
#biocLite("SPIA")
library(SPIA)
library("ROntoTools")
# kpg <- keggPathwayGraphs("hsa")
# kpg <- setEdgeWeights(kpg)
# kpg <- setNodeWeights(kpg, defaultWeight = 1)
# kpn <- keggPathwayNames("hsa")
# load("/wsu/home/gd/gd03/gd0393/Pathway/KEGG65.150Pathways.RData")
load("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/KEGG65.150Pathways.RData")
kpg <- setEdgeWeights(kpg)
kpg <- setNodeWeights(kpg, defaultWeight = 1)

getnodes=function(x){
  return (x@nodes);}


#generate the null distribution
#samePlat=c("GSE28146","GSE5281_EC","GSE5281_HIP","GSE5281_MTG","GSE5281_PC", "GSE5281_SFG", "GSE5281_VCX")
databases=c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
           "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")
#databases=c("GSE781", "GSE1297")
samePlat=as.list(databases)

runSPIA <- function (dataset) {
  dataset=samePlat[[2]]
  print(dataset)
  filename=paste(path,dataset,"/",dataset,"_annotation.txt",sep="")
  annotation=as.character((read.table(filename,header=F,sep="\t",stringsAsFactors=F))[1,1])
  keggnodes=lapply(kpg,getnodes)
  allGenes=unique(unlist(keggnodes))
  ano=paste(substring(annotation, 1, nchar(annotation)-3),"ENTREZID",sep="")
  require(annotation, character.only = TRUE)
  geneAno=as.data.frame(get(ano))
  geneAnno=geneAno[,2]
  names(geneAnno)=geneAno[,1]
  allProbes=unique(geneAno$probe_id)
  
  load(paste(path,dataset,"/",dataset,".RData",sep=""))      
  data <- get(paste("data_",dataset,sep=""))
  group <- get(paste("group_",dataset,sep=""))
  # data <- 2^data
  data <- data[rownames(data)%in%allProbes,]
  rownames(data) <- as.character(geneAnno[rownames(data)])
  data <- aggregate(data,by=list(rownames(data)),FUN=median)
  data <- data[paste("hsa:",data$Group.1,sep="")%in%allGenes,]
  rownames(data)=data$Group.1
  data <- data[,!colnames(data)%in%c("Group.1")]
  controlDat = data[,group$Group %in% 'c']
  diseaseDat = data[,group$Group %in% 'd']
  
  controlMean <- apply(controlDat,MARGIN=1,mean)
  diseaseMean <- apply(diseaseDat,MARGIN=1,mean)
  
  foldChange <- diseaseMean-controlMean
  names(foldChange) <- paste("hsa:",names(foldChange),sep="")
  
  pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
  
  names(pvalues) <- names(foldChange)
  DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
  DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
  DEGenes <- foldChange[names(DEGenes)]
  
  #peRes=pe(x=DEGenes,graphs=kpg,ref=names(foldChange),nboot=iterno)
  #temp <- Summary(peRes)
  #temp[,dim(temp)[2]+1] = rownames(temp)
  #rownames(temp)<-kpn[rownames(temp)]
  #colnames(temp)[dim(temp)[2]] = "pathway"
  #temp
  ########################
  names(DEGenes) <- gsub("hsa:", "", names(DEGenes))
  names(foldChange) <- gsub("hsa:", "", names(foldChange))
  
  res=spia(de=DEGenes, all=names(foldChange), organism="hsa",beta=NULL,nB=2000,plots=FALSE, verbose=TRUE,combine="fisher")
  
  res <- res[,-12]
  rownames(res) <- res[,1]
  res[,1] <- NULL
  res[,dim(res)[2]+1] <- res[,1]
  res[,dim(res)[2]] <- paste0("path:hsa",res[,dim(res)[2]])
  res[,1] <- NULL
  colnames(res)[dim(res)[2]] = "pathway"
  res
}

#library(multicore)
# resList <- mclapply(X=samePlat,FUN=runSPIA,mc.cores=5)
resList <- lapply(X=samePlat,FUN=runSPIA)

save(resList,file="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsSPIA.RData")

################ number of significant pathway ###########

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsSPIA.RData")
signPW = vector(mode = "numeric", length = 0)
for (i in 1:length(resList)) {
  signPW = c(signPW, length(which(resList[[i]]$pG<0.05))/150)
}

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/SPIA40Datasets.RData")
for (i in 1:length(myData)) {
  signPW = c(signPW, length(which(myData[[i]]$pG<0.05))/150)
}

mean(signPW)


