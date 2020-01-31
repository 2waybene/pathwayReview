library("GSEABase")
library("KEGG.db")
library("GOstats")

###############################################


setwd("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/");
path="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/"

filename="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/datasetslist.txt"
dataSets = as.character(t(read.table(filename,header=F,sep="\t",stringsAsFactors=F)))

PVAL <- 0.05
maxDE=400

load("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/KEGG65.15Pathways_Minh.RData") #load kpg, kpn, keggnodes
#load("/wsu/home/gd/gd03/gd0393/Pathway/KEGG65.150Pathways.RData")
load("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/GeneID2Symbol.RData") #load myGeneID, mapping ID to Symbols

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
loopNo=2000

generateNull = function (a) {
  it <- a[(nc+nd+1)]
  print(it)
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
    rownames(res) <- res$Term
    res <- res[,c("Term", "Pvalue")]
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

pathName <- kpn
resNullGOstats=data.frame(row.names=pathName, Name=pathName)
for (i in 1:loopNo) {
  # print(i)
  if (!is.na(resNull[[i]])) {
    if (length(resNull[[i]][,"Pvalue"])==0) { resNullGOstats[,paste("p_",i,sep="")] = NA }
    else {
      resNullGOstats[as.vector(resNull[[i]][,"Term"]),paste("p_",i,sep="")] = resNull[[i]][,"Pvalue"]
    }
  }
  
}

# save(resNullGOstats, file = "resNullGOstats.RData")

PDFPath = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/generateNullGOstats_color.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  temp = unlist(resNullGOstats[i,-1]) 
  temp = as.numeric(as.character(temp))
  hist(temp, breaks = 0:20/20, main = resNullGOstats[i,]$Name, col = "blue", xlab = "p-value")
}
dev.off()

file = paste0("GOstats",nc,"vs",nd,"_400_",2000,".RData",sep="") 
save(resNullGOstats,l,file=file)

################################# Pearson's Bias ################################# 
load(file = "GOstats15vs15_400_2000.RData")
library(sfsmisc)
path = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("pwNames.RData")

expected <- function(x) {
  y = density(x)
  # plot(y)
  expectedvalue = integrate.xy(y$x,y$x*y$y)
  expectedvalue
}

Pearson <- function(x) {
  # a = expected((x - expected(x))^3)/(sd(x)^3)
  a = expected((x-mean(x))^3)/(sd(x)^3)
  a
} 

bias = vector(mode="numeric", length=nrow(resNullGOstats))
names(bias) = resNullGOstats$Name

for (i in 1:nrow(resNullGOstats)) {
  x = as.numeric(resNullGOstats[i,-c(1:3)])
  x =  na.omit(x)
  bias[i] = Pearson(x)
}

bias0List = bias[which(bias >= 0.1)]
bias1List = bias[which(bias <= -0.1)]

bias0List = names(bias0List)
bias1List = names(bias1List)
bias0List %in% pwNames$Name
bias1List %in% pwNames$Name
length(bias0List) + length(bias1List)
nrow(resNullGOstats)

save(bias0List, bias1List, file = "GOstats_bias.RData")
