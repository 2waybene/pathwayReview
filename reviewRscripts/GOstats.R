library("GSEABase")
library("KEGG.db")
library("GOstats")


load("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/KEGG65.15Pathways_Minh.RData") #load kpg, kpn, keggnodes
pathName <- names(keggnodes)
mykeggframeData <- data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors=FALSE)
for (i in 1:length(keggnodes)) {
  a <- as.vector(keggnodes[[i]])
  for(j in 1:length(a)){
    mykeggframeData <- rbind(mykeggframeData, c(pathName[i], a[j]), stringsAsFactors=FALSE)
  }
}
colnames(mykeggframeData) <- c("Pathway", "gene")
mykeggframeData[,"Pathway"] <- gsub("path:hsa","",mykeggframeData[,"Pathway"])
mykeggframeData[,"gene"] <- gsub("hsa:","",mykeggframeData[,"gene"])

# mykeggFrame=KEGGFrame(mykeggframeData,organism="Homo sapiens")
# myGsc <- GeneSetCollection(mykeggFrame, setType = KEGGCollection())


filename="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/datasetslist.txt"
dataSets = as.character(t(read.table(filename,header=F,sep="\t",stringsAsFactors=F)))
filename="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/TargetPathwaysList.txt"
TargetPathways = as.character(t(read.table(filename,header=F,sep="\t", quote = "", stringsAsFactors=F)))

setwd("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/");
path="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/"


###############################################

PVAL <- 0.05
maxDE=400


GOstats <- function(dataset) {
  print(dataset)
  
  path="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/"
  load(paste(path,dataset,"/",dataset,".RData",sep=""))      
  data <- get(paste("gene_",dataset,sep=""))
  group <- get(paste("group_",dataset,sep=""))
  
  # data <- get(paste("gene_",dataset,sep=""))
  rownames(data) <- gsub("hsa:", "", rownames(data))
  data <- as.matrix(data)

  delRow <- sapply(seq(nrow(data)), function(x) { ifelse(sd(data[x,]) < 1e-10, FALSE, TRUE) })
  data <- data[delRow,]
  
  group <- get(paste("group_",dataset,sep=""))
  
  controlDat = data[,group$Group %in% 'c']
  diseaseDat = data[,group$Group %in% 'd']
  
  controlMean <- apply(controlDat,MARGIN=1,mean)
  diseaseMean <- apply(diseaseDat,MARGIN=1,mean)
  
  foldChange <- diseaseMean-controlMean
  
  pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)

  names(pvalues) <- names(foldChange)
  DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
  DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
  DEGenes <- foldChange[names(DEGenes)]
  
  genes <- names(DEGenes)
  
  universe <- names(pvalues)
  frame = mykeggframeData[mykeggframeData[,"gene"] %in% universe,]
  keggFrame = KEGGFrame(frame,organism="Homo sapiens")
  
  myGsc <- GeneSetCollection(keggFrame, setType = KEGGCollection())
  
  kparams <- GSEAKEGGHyperGParams(name="My Custom GSEA based annot Params", 
                                  geneSetCollection=myGsc, 
                                  geneIds = genes, 
                                  universeGeneIds = universe,  
                                  pvalueCutoff = 1, 
                                  testDirection = "over")
  kOver <- hyperGTest(kparams)
  res <- summary(kOver)
  res["Rank"] <- rank(res$Pvalue, ties.method="min");
  res
}
# dataSets = c("GSE15471", "GSE16515")
# subDataset <- dataSets[1:2]
result <- lapply(X = as.list(dataSets), FUN = GOstats)
# save(result, file = "WebGestaltResult.RData")

backingUp <- result
result <- backingUp

result <- lapply(X = result, function(x) { x <- x[,c("Term", "KEGGID", "Pvalue", "Rank")];} )

Rank <- vector("numeric")
pValue <- vector("numeric")
for (i in 1:length(TargetPathways)) {
  temp <- result[[i]]
  if (dim(temp[temp["Term"] == TargetPathways[i] ,])[1] == 0) {
    rank <- nrow(temp) + 1
    p <- NA
  } else {
    rank <- temp[temp["Term"] == TargetPathways[i] ,"Rank"]
    p <- temp[temp["Term"] == TargetPathways[i] ,"Pvalue"]
  }
  Rank <- c(Rank,rank)
  pValue <- c(pValue,p)
}

GOstats <- data.frame(Rank = Rank, pValue = pValue)
rownames(GOstats) <- dataSets
save(backingUp, result, GOstats, file = "GOstatsResult.RData")
