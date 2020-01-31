# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("EnrichmentBrowser", version = "3.8")
library("pROC")
library("EnrichmentBrowser")
library("WebGestaltR")
library(mGSZ)

path = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/mmu"
setwd(path)

#Create a gmt (gene matrix transposed) file
### Createss a gmt (gene matrix transposed) file such as those
### provided by mSigDB or geneSigDB, from an R list object.
### R list object that will be converted to GMT file.  Each element
### should contain a vector of gene names, and the names of the
### elements will used for the gene set names
myWriteGMT <- function (object, fname){
  if (class(object) != "list") stop("object should be of class 'list'")
  if(file.exists(fname)) unlink(fname)
  for (iElement in 1:length(object)){
    write.table(t(c(make.names(c(names(object)[iElement], NA)),object[[iElement]])),
                sep="\t",quote=FALSE,
                file=fname,append=TRUE,col.names=FALSE,row.names=FALSE)
  }
}

library("pROC")
library("org.Mm.eg.db")
xx <- as.list(org.Mm.egALIAS2EG)
for (i in 1:length(xx)) {
  if (xx[[i]] > 1) {
    xx[[i]] <- xx[[i]][1]
  }
}
mapping <- cbind(as.vector(names(xx)),as.vector(xx))
mapping <- apply(mapping, 2, toupper)
vectorMapping <- mapping[,1]
names(vectorMapping) <- mapping[,2]
myGeneID <- mapping[,1]
names(myGeneID) <- mapping[,2]
# load("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/mmu/kegg_mmu.RData")
# kpg = x$kpg
# kpn = x$kpn
# 
# getnodes=function(x){
#   return (x@nodes);}
# keggnodes=lapply(kpg,getnodes)
# save(kpg, kpn, keggnodes, file = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/mmu/KEGG65.15Pathways_mmu_Minh.RData")
load("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/mmu/KEGG65.15Pathways_mmu_Minh.RData")

pathNames <- kpn
names(pathNames) <- NULL


myKEGGgmt <- lapply(keggnodes, function(x) { x<-gsub("mmu:","",x); x<-myGeneID[x]; names(x)<-NULL; x<-na.omit(x); x<-as.vector(x);x})
referenceGenes <- unlist(myKEGGgmt, recursive = TRUE)
names(referenceGenes) <- NULL
referenceGenes <- unique(referenceGenes)
myWriteGMT(myKEGGgmt, "myKEGG_mmu.gmt")

###############################################

dataSets = c("GSE70302", "GSE70302", "GSE58120", "GSE46211", "GSE49166", "GSE50933", "GSE62999", "GSE57917", "GSE22873", "GSE6030", "GSE29048")
targetGenes = c("mmu:16175", "mmu:16176", "mmu:16428", "mmu:21813", "mmu:20893", "mmu:15903", "mmu:240672", "mmu:15379", "mmu:17874", "mmu:18012", "mmu:18609") 
targetSymbols = c("IL1a", "IL1b", "IL2", "TGFBR2", "BHLHE40", "ID3", "DUSP5", "ONECUR1", "MYD88", "NEUROD1", "PDX1")

# filename="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/datasetslist.txt"
# dataSets = as.character(t(read.table(filename,header=F,sep="\t",stringsAsFactors=F)))
# filename="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/TargetPathwaysList.txt"
# TargetPathways = as.character(t(read.table(filename,header=F,sep="\t", quote = "", stringsAsFactors=F)))

# path="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/mmu/WebGestalt"
# setwd(path);

geneFile <- "interestingGenes.txt"
refFile <- "referenceGenes.txt"
outputDirectory <- getwd()

# iterno=1000
PVAL <- 0.05
maxDE=400

myData <- list()
for (i in 1:length(dataSets)) {
  dataset = dataSets[i]
  targetID = targetGenes[i]
  targetSymbol = targetSymbols[i]
  
  print(dataset)
  filePath = paste0("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/mmu/dataset/", 
                    dataset,"_", targetSymbol, "/", dataset, "_", targetSymbol, ".RData")
  load(filePath)
    
  data <- get(paste("gene_",dataset,sep=""))
  group <- get(paste("group_",dataset,sep=""))
  
  rownames(data) <- gsub("mmu:", "", rownames(data))
  data <- as.matrix(data)
  rownames(data) <- myGeneID[rownames(data)]
  
  controlDat = data[,group$Group %in% 'c']
  diseaseDat = data[,group$Group %in% 'd']
  
  controlMean <- apply(controlDat,MARGIN=1,mean)
  diseaseMean <- apply(diseaseDat,MARGIN=1,mean)
  
  foldChange <- diseaseMean-controlMean
  #names(foldChange) <- paste("mmu:",names(foldChange),sep="")
  
  pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
    # for(j in 1:nrow(controlDat)) {
  #   print(j)
  #   t.test(t(controlDat[j,]),t(diseaseDat[j,]),paired=FALSE)
  # }
  
  names(pvalues) <- names(foldChange)
  DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
  DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
  DEGenes <- foldChange[names(DEGenes)]
  
  write(names(DEGenes), file = "interestingGenes.txt")
  write(intersect(referenceGenes, rownames(data)) , file = "referenceGenes.txt")
  
  
  enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                              enrichDatabase="others", enrichDatabaseFile = "myKEGG_mmu.gmt",
                              interestGeneFile=geneFile, enrichDatabaseType = "genesymbol",
                              fdrThr = 1, topThr = 150, reportNum = 150, sigMethod = "top",
                              setCoverNum = 150,minNum = 1, maxNum = 4000,
                              interestGeneType="genesymbol", referenceGeneFile=refFile,
                              referenceGeneType="genesymbol", isOutput=TRUE,
                              outputDirectory=outputDirectory, projectName=paste0("WebGestaltORA",dataset,targetSymbol))
  
  res <- as.matrix(read.table(file = paste0("Project_WebGestaltORA",dataset,targetSymbol,
                                                 "/enrichment_results_WebGestaltORA",dataset,targetSymbol,".txt"), sep="\t", header=TRUE))
  myData[[i]] <- res
}

save(myData, file = "WebGestaltResult_mmu.RData")

backingUp <- myData
myData <- backingUp

myData <- lapply(X = myData, function(x) { 
  x <- as.data.frame(x);
  x <- x[,c("geneSet", "pValue")];
  x[,1] <- gsub("path.", "path:", x[,1]); 
  x$pValue <- as.numeric(as.character(x$pValue));
  x["PathwayName"] <- kpn[x[,1]];
  x["Rank"] <- rank(x$pValue,ties.method="min");
  x;
} )



TP <- vector("numeric")
FP <- vector("numeric")
FN <- vector("numeric")
TN <- vector("numeric")
sensitivity <- vector("numeric")
specificity <- vector("numeric")
accuracy <- vector("numeric")
auc <- vector("numeric")

for (i in 1:length(myData)) {
  print(i)
  targetID <- targetGenes[i]
  mylist <- lapply(keggnodes, function(x) { length(intersect(x, targetID))} )
  mylist = as.vector(mylist)
  names(mylist) = names(keggnodes)
  posPredict = mylist[which(mylist == 1)]
  negPredict = mylist[!mylist %in% posPredict]
  
  temp <- myData[[i]]
  # pwID <- as.vector(temp["geneSet"])
  # pwName <- as.vector(temp[,"PathwayName"])
  # pvalue <- as.vector(temp[,"pValue"])
  # temp = cbind(pwName,pvalue)
  # rownames(temp) <- pwID

  # r <- nrow(temp)
  # for (j in 1: length(posPredict)){
  #   rn <- rownames(temp)
  #   temp <- rbind(temp, unlist(c(names(posPredict[j]), 1, unname(kpn[names(posPredict[j])]), r+1)))
  #   rn <- c(rn, names(posPredict[j]))
  # }
  label <- rep(NA, nrow(temp))
  names(label) <- temp[,"geneSet"]
  label[names(label)%in%names(posPredict)] = 1
  label[names(label)%in%names(negPredict)] = 0
  
  temp = cbind.data.frame(temp,label)
  colnames(temp) = c("geneSet", "pValue", "PathwayName", "Rank", "label")
  temp = na.omit(temp)
  temp <- temp[order(temp$pValue),]
  labels = factor(label)
  
  if(nlevels(labels) > 1) {
    cutoff <- 0.05
    TP[i] = nrow(temp[as.numeric(temp[,"pValue"]) <= cutoff & temp[,"label"] == 1, ])
    FP[i] = nrow(temp[as.numeric(temp[,"pValue"]) <= cutoff & temp[,"label"] == 0, ])
    FN[i] = nrow(temp[as.numeric(temp[,"pValue"]) > cutoff & temp[,"label"] == 1, ])
    TN[i] = nrow(temp[as.numeric(temp[,"pValue"]) > cutoff & temp[,"label"] == 0, ])
    
    # Calculate sensitivity
    sensitivity[i] = TP[i]/(TP[i] + FN[i])
    sensitivity[i]
    accuracy[i] = (TN[i] + TP[i])/(TP[i] + FP[i] + FN[i] + TN[i])
    accuracy[i]
    specificity[i] = TN[i]/(TN[i]+FP[i])
    specificity[i]
    
    rocobj <- roc(labels, as.numeric(temp[,"pValue"]),percent = T)  
    plot(rocobj,legacy.axes=TRUE, main= "ROC curve of WebGestalt", percent=TRUE, col="red", print.auc=T )
    auc(roc(labels, as.numeric(temp[,"pValue"])))[1]
    auc[i] <- auc(roc(labels, as.numeric(temp[,"pValue"])))[1]
  } else {
    auc[i] <- 0.5
  }
  
}

save(auc, sensitivity, accuracy, specificity, file = "WebGestaltResult_mmu.RData")








