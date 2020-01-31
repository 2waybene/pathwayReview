if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnrichmentBrowser", version = "3.8")

library("EnrichmentBrowser")
library("WebGestaltR")
library(mGSZ)


###############################################


# setwd("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/WebGestaltNull");
# path="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/WebGestaltNull"
setwd("/wsu/home/gd/gd03/gd0393/Pathway/WebGestaltNull");
path="/wsu/home/gd/gd03/gd0393/Pathway/WebGestaltNull"


filename="/wsu/home/gd/gd03/gd0393/Pathway/dataset/datasetslist.txt"
dataSets = as.character(t(read.table(filename,header=F,sep="\t",stringsAsFactors=F)))
filename="/wsu/home/gd/gd03/gd0393/Pathway/dataset/TargetPathwaysList.txt"
TargetPathways = as.character(t(read.table(filename,header=F,sep="\t", quote = "", stringsAsFactors=F)))


geneFile <- "interestingGenes.txt"
refFile <- "referenceGenes.txt"
outputDirectory <- getwd()

iterno=1000
PVAL <- 0.05
maxDE=400

# load("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/KEGG65.15Pathways_Minh.RData") #load kpg, kpn, keggnodes
load("/wsu/home/gd/gd03/gd0393/Pathway/KEGG65.15Pathways_Minh.RData")
load("/wsu/home/gd/gd03/gd0393/Pathway/GeneID2Symbol.RData") #load myGeneID, mapping ID to Symbols

pathNames <- kpn
names(pathNames) <- NULL
myKEGGgmt <- lapply(keggnodes, function(x) { x<-gsub("hsa:","",x); x<-myGeneID[x]; names(x)<-NULL; x<-na.omit(x); x<-as.vector(x);x})
referenceGenes <- unlist(myKEGGgmt, recursive = TRUE)
names(referenceGenes) <- NULL
referenceGenes <- unique(referenceGenes)

###################### generate the null distribution ###################### 

# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/mygeneSym.RData")
#load("/wsu/home/gd/gd03/gd0393/Pathway/mygeneSym.RData")

Samples <- NULL
intersectGenes <- NULL
for (i in 1:length(dataSets)) {
  dataset=dataSets[i]
  print(c("i: ", i, "; dataset: ", dataset))
  path="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/"
  load(paste(path,dataset,"/",dataset,".RData",sep=""))  
  data <- get(paste("gene_",dataset,sep=""))
  group <- get(paste("group_",dataset,sep=""))
  cn = group[group$Group=="c","Sample"]
  data <- data[,cn]
  
  rownames(data) <- gsub("hsa:", "", rownames(data))
  data <- as.matrix(data)
  rownames(data) <- myGeneID[rownames(data)]
  data <- aggregate(data,by=list(rownames(data)),FUN=median)
  rownames(data)=data$Group.1
  data <- data[,-1]  
  
  if (is.null(intersectGenes)) {
    intersectGenes = rownames(data)
  } else {
    intersectGenes = intersect(intersectGenes, rownames(data))
  }
  if (is.null(Samples)) {
    Samples = data
  } else {
    Samples <- cbind(Samples[intersectGenes,],data[intersectGenes,])
  }
  
}

nc=15
nd=15
#loopNo=5000
loopNo=3

generateNull = function (a) {
  it <- a[(nc+nd+1)]
  a <- a[-(nc+nd+1)]
  controlDat = Samples[,a[1:nc]]
  diseaseDat = Samples[,a[(nc+1):(nc+nd)]]
  #wholeDat <- cbind(controlDat, diseaseDat)
  ##########
  
  controlMean <- apply(controlDat,MARGIN=1,mean)
  diseaseMean <- apply(diseaseDat,MARGIN=1,mean)
  
  foldChange <- diseaseMean-controlMean
  # names(foldChange) <- paste("hsa:",names(foldChange),sep="")
  
  pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
  
  names(pvalues) <- names(foldChange)
  
  DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
  if (length(DEGenes) != 0) {
    DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
    DEGenes <- foldChange[names(DEGenes)]
    write(names(DEGenes), file = "interestingGenes.txt")
    write(intersect(referenceGenes, rownames(data)) , file = "referenceGenes.txt")
    
    if (all(names(DEGenes) %in% intersect(referenceGenes, rownames(data)) == TRUE)) {
      enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                  enrichDatabase="others", enrichDatabaseFile = "myKEGG.gmt",
                                  interestGeneFile=geneFile, enrichDatabaseType = "genesymbol",
                                  fdrThr = 1, topThr = 150, reportNum = 150, sigMethod = "top",
                                  setCoverNum = 150,minNum = 1, maxNum = 4000,
                                  interestGeneType="genesymbol", referenceGeneFile=refFile,
                                  referenceGeneType="genesymbol", isOutput=TRUE,
                                  outputDirectory=outputDirectory, projectName=paste0("WebGestaltORA", it))
      
      res <- as.matrix(read.table(file = paste0("Project_WebGestaltORA",it,
                                                "/enrichment_results_WebGestaltORA",it,".txt"), sep="\t", header=TRUE))
      res <- as.data.frame(res);
      res <- res[,c("geneSet", "pValue")];
      res[,1] <- gsub("path.", "path:", res[,1]); 
      res$pValue <- as.numeric(as.character(res$pValue));
      res["PathwayName"] <- kpn[res[,1]];
      res
    } else { res <- NA }
    
  } else { res <- NA  }
  res
}

l <- list()
for (i in 1:loopNo) {
  a=sample(colnames(Samples),nc+nd,replace=FALSE)
  a <- c(a,i)
  l[[i]] <- a
}


resNull <- lapply(X=l,FUN=generateNull)

resNullWebGestalt=data.frame(row.names=pathNames, Name=pathNames)
for (i in 1:loopNo) {
  if (length(resNull[[i]][,"pValue"])==0) { resNullWebGestalt[,paste("p_",i,sep="")] = NA }
  else {
    resNullWebGestalt[as.vector(resNull[[i]][,"PathwayName"]),paste("p_",i,sep="")] = resNull[[i]][,"pValue"]
  }
}

save(resNullWebGestalt, file = "resNullWebGestalt.RData")
