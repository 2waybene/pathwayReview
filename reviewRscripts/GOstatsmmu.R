library("GSEABase")
library("KEGG.db")
library("GOstats")
library("pROC")

path = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/mmu"
setwd(path)

load("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/mmu/KEGG65.15Pathways_mmu_Minh.RData")
pathName <- names(keggnodes)
names(kpn) <- gsub("path:mmu", "", names(kpn))
mykeggframeData <- data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors=FALSE)
for (i in 1:length(keggnodes)) {
  a <- as.vector(keggnodes[[i]])
  for(j in 1:length(a)){
    mykeggframeData <- rbind(mykeggframeData, c(pathName[i], a[j]), stringsAsFactors=FALSE)
  }
}
colnames(mykeggframeData) <- c("Pathway", "gene")
mykeggframeData[,"Pathway"] <- gsub("path:mmu","",mykeggframeData[,"Pathway"])
mykeggframeData[,"gene"] <- gsub("mmu:","",mykeggframeData[,"gene"])


###############################################

dataSets = c("GSE70302", "GSE70302", "GSE58120", "GSE46211", "GSE49166", "GSE50933", "GSE62999", "GSE57917", "GSE22873", "GSE6030", "GSE29048")
targetGenes = c("mmu:16175", "mmu:16176", "mmu:16428", "mmu:21813", "mmu:20893", "mmu:15903", "mmu:240672", "mmu:15379", "mmu:17874", "mmu:18012", "mmu:18609") 
targetSymbols = c("IL1a", "IL1b", "IL2", "TGFBR2", "BHLHE40", "ID3", "DUSP5", "ONECUR1", "MYD88", "NEUROD1", "PDX1")

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
  delRow <- sapply(seq(nrow(data)), function(x) { ifelse(sd(data[x,]) < 1e-10, FALSE, TRUE) })
  data <- data[delRow,]
  
  controlDat = data[,group$Group %in% 'c']
  diseaseDat = data[,group$Group %in% 'd']
  
  controlMean <- apply(controlDat,MARGIN=1,mean)
  diseaseMean <- apply(diseaseDat,MARGIN=1,mean)
  
  foldChange <- diseaseMean-controlMean
  #names(foldChange) <- paste("mmu:",names(foldChange),sep="")
  
  pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
  names(pvalues) <- names(foldChange)
  DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
  DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
  DEGenes <- foldChange[names(DEGenes)]
  
  genes <- names(DEGenes)
  
  universe <- names(pvalues)
  frame = mykeggframeData[mykeggframeData[,"gene"] %in% universe,]
  keggFrame = KEGGFrame(frame,organism="Mus musculus")
  
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
  res[,"Term"] <- kpn[res[,"KEGGID"]]
  myData[[i]] <- res
}

save(myData, file = "GOstatsResult_mmu.RData")

backingUp <- myData
myData <- backingUp

myData <- lapply(X = myData, function(x) { 
  rownames(x) <- x[,"KEGGID"]
  x <- x[,c("Term", "Pvalue", "Rank")];
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
  names(mylist) <- gsub("path:mmu", "", names(mylist))
  posPredict = mylist[which(mylist == 1)]
  negPredict = mylist[!mylist %in% posPredict]
  
  temp <- myData[[i]]

  label <- rep(NA, nrow(temp))
  names(label) <- rownames(temp)
  label[names(label)%in%names(posPredict)] = 1
  label[names(label)%in%names(negPredict)] = 0
  
  temp = cbind.data.frame(temp,label)
  colnames(temp) = c("Term", "Pvalue", "Rank", "label")
  temp = na.omit(temp)
  # temp <- temp[order(temp$Pvalue),]
  labels = factor(label)
  
  if(nlevels(labels) > 1) {
    cutoff <- 0.05
    TP[i] = nrow(temp[as.numeric(temp[,"Pvalue"]) <= cutoff & temp[,"label"] == 1, ])
    FP[i] = nrow(temp[as.numeric(temp[,"Pvalue"]) <= cutoff & temp[,"label"] == 0, ])
    FN[i] = nrow(temp[as.numeric(temp[,"Pvalue"]) > cutoff & temp[,"label"] == 1, ])
    TN[i] = nrow(temp[as.numeric(temp[,"Pvalue"]) > cutoff & temp[,"label"] == 0, ])
    
    # Calculate sensitivity
    sensitivity[i] = TP[i]/(TP[i] + FN[i])
    sensitivity[i]
    accuracy[i] = (TN[i] + TP[i])/(TP[i] + FP[i] + FN[i] + TN[i])
    accuracy[i]
    specificity[i] = TN[i]/(TN[i]+FP[i])
    specificity[i]
    
    rocobj <- roc(labels, as.numeric(temp[,"Pvalue"]),percent = T)  
    # plot(rocobj,legacy.axes=TRUE, main= "ROC curve of GOstats", percent=TRUE, col="red", print.auc=T )
    auc(roc(labels, as.numeric(temp[,"Pvalue"])))[1]
    auc[i] <- auc(roc(labels, as.numeric(temp[,"Pvalue"])))[1]
  } else {
    auc[i] <- 0.5
  }
  
}

save(auc, sensitivity, accuracy, specificity, file = "GOstatsResult_mmu.RData")

