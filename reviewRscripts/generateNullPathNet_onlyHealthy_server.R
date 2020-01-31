setwd("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/")
path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/"
# setwd("/wsu/home/gd/gd03/gd0393/Pathway/dataset/");
# path="/wsu/home/gd/gd03/gd0393/Pathway/dataset/"
 source("http://bioconductor.org/biocLite.R")
 biocLite("PathNet")
 biocLite("PathNetData")
iterno=2000
PVAL <- 0.05
maxDE=400

#library("ROntoTools")
library("PathNet")
library("PathNetData")

#load("/wsu/home/gd/gd03/gd0393/Pathway/KEGG65.150Pathways.RData")
#kpg <- setEdgeWeights(kpg)
#kpg <- setNodeWeights(kpg, defaultWeight = 1)

current <- getwd()
setwd(system.file(dir="extdata", package="PathNetData"))
brain_regions <- as.matrix(read.table(file = "brain_regions_data.txt", sep = "\t", header = T))
disease_progression <- as.matrix(read.table(file = "disease_progression_data.txt", sep = "\t", header = T))
A <- as.matrix(read.table(file = "adjacency_data.txt", sep = "\t", header = T))
pathway <- read.table(file = "pathway_data.txt", sep = "\t", header = T)
setwd(current)

#getnodes=function(x){
#  return (x@nodes);}


#generate the null distribution
samePlat=c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
          "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
          "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
          "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
          "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
          "GSE24739_G1","GSE32676","GSE4183","GSE7305")
#samePlat=c("GSE5281_EC", "GSE5281_HIP","GSE5281_VCX")
Samples <- NULL
intersectGenes <- NULL
for (i in 1:length(samePlat)) {
  dataset = samePlat[i]
  print(dataset)
  load(paste(path, dataset,"/",dataset,".RData",sep=""))
  
  annotation=get(paste("annotation_",dataset,sep=""))
  #design=get(paste("design_",dataset,sep=""))
  design=c("Not Paired")
  
  ano=paste(substring(annotation, 1, nchar(annotation)-3),"ENTREZID",sep="")
  require(annotation, character.only = TRUE)
  geneAno=as.data.frame(get(ano))
  geneAnno=geneAno[,2]
  names(geneAnno)=geneAno[,1]
  allProbes=unique(geneAno$probe_id)
  
  data <- get(paste("data_",dataset,sep=""))
  group <- get(paste("group_",dataset,sep=""))
  
  
  
  data <- 2^get(paste("data_",dataset,sep=""))
  data <- data[rownames(data)%in%allProbes,]
  rownames(data) <- as.character(geneAnno[rownames(data)])
  data <- aggregate(data,by=list(rownames(data)),FUN=median)
  #mydat <- mydat[paste("hsa:",mydat$Group.1,sep="")%in%allGenes,]
  
  cn = rownames(group[group$Group=="c",])
  data <- data[,cn]
  
  
  
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
loopNo=2000

generateNull = function (a) {
  controlDat = Samples[,a[1:nc]]
  
  diseaseDat = Samples[,a[(nc+1):(nc+nd)]]
  #wholeDat <- cbind(controlDat, diseaseDat)
  ##########
  
  controlMean <- log(apply(controlDat,MARGIN=1,mean),base=2)
  diseaseMean <- log(apply(diseaseDat,MARGIN=1,mean),base=2)
  
  foldChange <- diseaseMean-controlMean
  #names(foldChange) <- paste("hsa:",names(foldChange),sep="")
  
  pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
  
  names(pvalues) <- names(foldChange)
  
  pvalueMatrix <- matrix(c(as.numeric(names(pvalues)),as.numeric(pvalues)), ncol = 2)
  colnames(pvalueMatrix) <- c("Gene.ID", "EC")
  pvalueMatrix[,2] <- -log10(pvalueMatrix[,2])
  
  
  gene_ID <- pvalueMatrix[,1]
  
  A <- A[rownames(A) %in% gene_ID, rownames(A) %in% gene_ID]
  
  
  
  results <- PathNet(Enrichment_Analysis = TRUE,
                     DirectEvidence_info = pvalueMatrix[1:12000,],
                     Adjacency = A,
                     pathway = pathway,
                     Column_DirectEvidence = 2,
                     n_perm = 2000, threshold = 0.05)
  temp <- results$enrichment_results
  temp
}

l <- list()
for (i in 1:loopNo) {
  a=sample(colnames(Samples),nc+nd,replace=FALSE)
  l[[i]] <- a
}

#library(multicore)
resNull <- mclapply(X=l,FUN=generateNull,mc.cores=10)
#resNull <- lapply(X=l,generateNull)

pathways <- resNull[[1]]$Name

#save(resNull, file="/Users/GaMinh/PathNet.RData")
#load("/Users/GaMinh/PathNet.RData")
resNullPathNet=data.frame(row.names=pathways, Name=pathways,M=rep(NA,length(pathways)),std=rep(NA,length(pathways)))
for (i in 1:loopNo) {
  if (length(resNull[[i]]$p_PathNet)==0) {resNullPathNet[,paste("p_PathNet",i,sep="")]=NA}
  else {resNullPathNet[as.vector(resNull[[i]]$Name),paste0("p_PathNet",i)]=resNull[[i]]$p_PathNet}
}
resNullPathNet$M=apply(resNullPathNet[4:(loopNo+3)],1,FUN=mean,na.rm=T)
resNullPathNet$std=apply(resNullPathNet[4:(loopNo+3)],1,FUN=sd,na.rm=T)

file = paste("PathNet",nc,"vs",nd, "_",loopNo,".RData",sep="") 
save(resNullPathNet,l,file=file)

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/PathNet15vs15_2000(2nd run).RData")

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullPathNet.pdf"
pathName <- resNullPathNet[,"Name"]
pdf(file=PDFPath)

for (i in seq(length(pathName)))  {

    print(i)
    hist(as.numeric(resNullPathNet[i,-c(1,2,3)]), breaks = 50, main = resNullPathNet[i,]$Name)
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullPathNet_color.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
    print(i)
    hist(as.numeric(resNullPathNet[i,-c(1,2,3)]), breaks = 50, main = resNullPathNet[i,]$Name, col = "blue")
}
dev.off()

