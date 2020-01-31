# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("EnrichmentBrowser", version = "3.8")

library("EnrichmentBrowser")
library("WebGestaltR")
library(mGSZ)

#Create a gmt (gene matrix transposed) file
### Createss a gmt (gene matrix transposed) file such as those
### provided by mSigDB or geneSigDB, from an R list object.
### R list object that will be converted to GMT file.  Each element
### should contain a vector of gene names, and the names of the
### elements will used for the gene set names
# myWriteGMT <- function (object, fname){
#   if (class(object) != "list") stop("object should be of class 'list'")
#   if(file.exists(fname)) unlink(fname)
#   for (iElement in 1:length(object)){
#     write.table(t(c(make.names(c(names(object)[iElement], NA)),object[[iElement]])),
#                 sep="\t",quote=FALSE,
#                 file=fname,append=TRUE,col.names=FALSE,row.names=FALSE)
#   }
# }

############# mapping ##################################
library("org.Hs.eg.db")

# mapping gene ID to gene symbol
x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
mapping <- cbind(as.vector(names(xx)),as.vector(xx))
mapping <- apply(mapping, 2, toupper)
colnames(mapping) <- c("gene_annotation", "gene_id")
myGeneID = unlist(mapping[,"gene_id"])
names(myGeneID) = unlist(mapping[,"gene_annotation"])

###############################################

filename="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/datasetslist.txt"
dataSets = as.character(t(read.table(filename,header=F,sep="\t",stringsAsFactors=F)))
filename="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/TargetPathwaysList.txt"
TargetPathways = as.character(t(read.table(filename,header=F,sep="\t", quote = "", stringsAsFactors=F)))

setwd("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/WebGestalt");
path="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/WebGestalt"


geneFile <- "interestingGenes.txt"
refFile <- "referenceGenes.txt"
outputDirectory <- getwd()

iterno=1000
PVAL <- 0.05
maxDE=400

load("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/KEGG65.15Pathways_Minh.RData") #load kpg, kpn, keggnodes

pathNames <- kpn
names(pathNames) <- NULL


myKEGGgmt <- lapply(keggnodes, function(x) { x<-gsub("hsa:","",x); x<-myGeneID[x]; names(x)<-NULL; x<-na.omit(x); x<-as.vector(x);x})
referenceGenes <- unlist(myKEGGgmt, recursive = TRUE)
names(referenceGenes) <- NULL
referenceGenes <- unique(referenceGenes)
# myWriteGMT(myKEGGgmt, "myKEGG.gmt")

WebGestalt <- function(dataset) {
  print(dataset)
  
  path="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/"
  load(paste(path,dataset,"/",dataset,".RData",sep=""))      
  data <- get(paste("gene_",dataset,sep=""))
  group <- get(paste("group_",dataset,sep=""))
  
  # data <- get(paste("gene_",dataset,sep=""))
  rownames(data) <- gsub("hsa:", "", rownames(data))
  data <- as.matrix(data)
  rownames(data) <- myGeneID[rownames(data)]
  data <- aggregate(data,by=list(rownames(data)),FUN=median)
  rownames(data)=data$Group.1
  data <- data[,-1]
  delRow <- sapply(seq(nrow(data)), function(x) { ifelse(sd(data[x,]) < 1e-10, FALSE, TRUE) })
  data <- data[delRow,]
  
  group <- get(paste("group_",dataset,sep=""))
  
  controlDat = data[,group$Group %in% 'c']
  diseaseDat = data[,group$Group %in% 'd']
  
  controlMean <- apply(controlDat,MARGIN=1,mean)
  diseaseMean <- apply(diseaseDat,MARGIN=1,mean)
  
  foldChange <- diseaseMean-controlMean
  
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
                              enrichDatabase="others", enrichDatabaseFile = "myKEGG.gmt",
                              interestGeneFile=geneFile, enrichDatabaseType = "genesymbol",
                              fdrThr = 1, topThr = 150, reportNum = 150, sigMethod = "top",
                              setCoverNum = 150,minNum = 1, maxNum = 4000,
                              interestGeneType="genesymbol", referenceGeneFile=refFile,
                              referenceGeneType="genesymbol", isOutput=TRUE,
                              outputDirectory=outputDirectory, projectName=paste0("WebGestaltORA",dataset))
  
  res <- as.matrix(read.table(file = paste0("Project_WebGestaltORA",dataset,
                                                 "/enrichment_results_WebGestaltORA",dataset,".txt"), sep="\t", header=TRUE))
  res
}
# dataSets = c("GSE15471", "GSE16515")
result <- lapply(X = as.list(dataSets), FUN = WebGestalt)
# save(result, file = "WebGestaltResult.RData")

backingUp <- result
result <- backingUp

result <- lapply(X = result, function(x) { 
  x <- as.data.frame(x);
  x <- x[,c("geneSet", "pValue")];
  x[,1] <- gsub("path.", "path:", x[,1]); 
  x$pValue <- as.numeric(as.character(x$pValue));
  x["PathwayName"] <- kpn[x[,1]];
  x["Rank"] <- rank(x$pValue,ties.method="min");
  x;
  } )

Rank <- vector("numeric")
pValue <- vector("numeric")
for (i in 1:length(TargetPathways)) {
  temp <- result[[i]]
  if (dim(temp[temp["PathwayName"] == TargetPathways[i] ,])[1] == 0) {
    rank <- nrow(temp) + 1
    p <- NA
  } else {
    rank <- temp[temp["PathwayName"] == TargetPathways[i] ,"Rank"]
    p <- temp[temp["PathwayName"] == TargetPathways[i] ,"pValue"]
  }
  Rank <- c(Rank,rank)
  pValue <- c(pValue,p)
}

WebGestaltResult <- data.frame(Rank = Rank, pValue = pValue)
rownames(WebGestaltResult) <- dataSets
save(backingUp, WebGestaltResult, file = "WebGestaltResult.RData")
