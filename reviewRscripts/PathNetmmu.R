# library(limma)
# y = getGeneKEGGLinks(species.KEGG = "mmu", convert = FALSE)
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/HighEdgeS/inst/extdata/kpgNames_mmu_UpdatedPathways.RData")
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/data/kpg_mmu_UpdatedPathways.RData")
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/data/KEGGgraph_mmu_UpdatedPathways.RData")
load("/Users/GaMinh/Downloads/kegg_mmu.RData")
load("/Users/GaMinh/Downloads/altEntrez_GSE22873_Myd88_mmu_SAMER.RData")
kpg = x$kpg
kpn = x$kpn

getnodes=function(x){
  return (x@nodes);}
keggnodes=lapply(kpg,getnodes)

# source("https://bioconductor.org/biocLite.R")
# biocLite("BLMA")
# library("BLMA")
library("ROntoTools")
# x <- loadKEGGPathways(organism = "mmu", updateCache = T)
iterno=2000
PVAL <- 0.05
maxDE=400
mydata = DataObject$exprTable
rownames(mydata) = paste0("mmu:",rownames(mydata))
foldChange = mydata$logFC
names(foldChange) = rownames(mydata)
pvalues = mydata$P.Value
names(pvalues) = rownames(mydata)

DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)

DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
DEGenes <- foldChange[names(DEGenes)]

pvalueMatrix <- matrix(c(as.numeric(names(pvalues)),as.numeric(pvalues)), ncol = 2)
colnames(pvalueMatrix) <- c("Gene.ID", "EC")
pvalueMatrix[,2] <- -log10(pvalueMatrix[,2])




library("PathNet")
library("PathNetData")

# load the matrix A for pathnet
setwd(system.file(dir="extdata", package="PathNetData"))
A <- as.matrix(read.table(file = "adjacency_data.txt", sep = "\t", header = T))
pathway <- read.table(file = "pathway_data.txt", sep = "\t", header = T)

# Set working directory
path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)

dataSets=c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
           "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")

myData <- list()
#for (i in 1:length(dataSets)) {
for (i in 1:6) {
  dataset = dataSets[i]
  load(paste(path,"dataset/", dataset,"/",dataset,".RData",sep=""))
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
  mydat <- aggregate(data,by=list(rownames(data)),FUN=min)
  #mydat <- aggregate(data,by=list(rownames(data)),FUN=median)
  #mydat <- mydat[paste("hsa:",mydat$Group.1,sep="")%in%allGenes,]
  
  controlDat = mydat[,rownames(group[group$Group=="c",])]
  rownames(controlDat) = as.character(mydat$Group.1)
  
  diseaseDat = mydat[,rownames(group[group$Group=="d",])]
  rownames(diseaseDat) = as.character(mydat$Group.1)
  
  
  controlMean <- log(apply(controlDat,MARGIN=1,mean),base=2)
  diseaseMean <- log(apply(diseaseDat,MARGIN=1,mean),base=2)
  
  foldChange <- diseaseMean-controlMean
  
  pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
  names(pvalues) <- names(foldChange)
  
  
  pvalueMatrix <- matrix(c(as.numeric(names(pvalues)),as.numeric(pvalues)), ncol = 2)
  colnames(pvalueMatrix) <- c("Gene.ID", "EC")
  pvalueMatrix[,2] <- -log10(pvalueMatrix[,2])
  
  
  gene_ID <- pvalueMatrix[,1]
  
  A <- A[rownames(A) %in% gene_ID, rownames(A) %in% gene_ID]
  
  results <- PathNet(Enrichment_Analysis = TRUE,
                     DirectEvidence_info = pvalueMatrix[1:2000,],
                     Adjacency = A,
                     pathway = pathway,
                     Column_DirectEvidence = 2,
                     n_perm = 2000, threshold = 0.05)
  
  myData[[i]] <- results$enrichment_results
}
names(myData) <- dataSets
save(myData,file=paste(path,"alldatasetsPathNet.RData",sep=""))