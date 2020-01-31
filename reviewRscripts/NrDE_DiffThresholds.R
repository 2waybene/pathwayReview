library("pROC")
library("org.Mm.eg.db")
load("/Users/GaMinh/Downloads/kegg_mmu.RData")
kpg = x$kpg
kpn = x$kpn

getnodes=function(x){
  return (x@nodes);}
keggnodes=lapply(kpg,getnodes)
iterno=2000
PVAL <- 0.05
maxDE=400

path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/mmu/"
setwd(path);
load("/Users/GaMinh/Downloads/kegg_mmu.RData")

datasets = c("GSE70302", "GSE70302", "GSE58120", "GSE46211", "GSE49166", "GSE50933", "GSE62999", "GSE57917")
targetGenes = c("mmu:16175", "mmu:16176", "mmu:16428", "mmu:21813", "mmu:20893", "mmu:15903", "mmu:240672", "mmu:15379") 
targetGenesSymbol = c("IL1a", "IL1b", "IL2", "TGFBR2", "BHLHE40", "ID3", "DUSP5", "ONECUR1")

methodDE = c("Classic", "400DE")
j = 1 #classic or 400DE?
thresholdFC = 0.5
nrDE = rep(NA, length(datasets))

for (i in 1:length(datasets)) {
  dataset = datasets[i]
  targetID = targetGenes[i]
  targetSymbol = targetGenesSymbol[i]
  filePath = paste0(path, "dataset/" , dataset,"_", targetSymbol, "/", dataset, "_", targetSymbol, ".RData")
  load(filePath)
  
  data <- get(paste("gene_",dataset,sep=""))
  group <- get(paste("group_",dataset,sep=""))
  
  controlDat = data[,group$Group %in% 'c']
  diseaseDat = data[,group$Group %in% 'd']
  
  controlMean <- apply(controlDat,MARGIN=1,mean)
  diseaseMean <- apply(diseaseDat,MARGIN=1,mean)
  
  foldChange <- diseaseMean-controlMean
  #names(foldChange) <- paste("mmu:",names(foldChange),sep="")
  
  pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
  names(pvalues) <- names(foldChange)
  
  if (j == 2) {
    ######## 400 DE GENES METHOD ######## 
    print("400DE")
    DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
    DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
    DEGenes <- foldChange[names(DEGenes)]
  } else {
    ######## TRADITIONAL WAY TO CHOOSE DE GENES ########
    print("Classic")
    adj.pvalues <- p.adjust(pvalues, method = 'fdr')
    DEGenes <- foldChange[which(abs(foldChange[adj.pvalues<PVAL]) > thresholdFC)]
  }
  
  nrDE[i] = length(DEGenes)
}
nrDE



################ SAMER's DATASETS #################
path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/mmu/"
setwd(path);
load("/Users/GaMinh/Downloads/kegg_mmu.RData")

filePath = c("/Users/GaMinh/Downloads/altEntrez_GSE22873_Myd88_mmu_SAMER.RData", 
             "Code and Data/data/altEntrez_GSE6030_Neurod1_mmu_SAMER.RData",
             "Code and Data/data/altEntrez_GSE29048_Pdx1_mmu_SAMER.RData")
EsetFile = c("Code and Data/data/DatasetsGAGE/GSE22873_NORM/GSE22873_altCdf.RData",
             "Code and Data/data/DatasetsGAGE/GSE6030_NORM/GSE6030_altCdf.RData",
             "Code and Data/data/DatasetsGAGE/GSE29048_NORM/GSE29048_altCdf.RData")

datasets = c("GSE22873_Myd88","GSE6030NEUROD1", "GSE29048_Pdx1")
targetGene = c("mmu:17874", "mmu:18012", "mmu:18609") 
methodDE = c("Classic", "400DE")
j = 1
thresholdFC = 0.5

for (i in 1:length(datasets)) {
  load(filePath[i])
  load(EsetFile[i])
  dataset = datasets[i]
  target = targetGene[i]
  
  iterno=2000
  PVAL <- 0.05
  maxDE=400
  mydata = DataObject$exprTable
  rownames(mydata) = paste0("mmu:",rownames(mydata))
  foldChange = mydata$logFC
  names(foldChange) = rownames(mydata)
  pvalues = mydata$P.Value
  names(pvalues) = rownames(mydata)
  
  if (j == 2) {
    ######## 400 DE GENES METHOD ######## 
    #print("400DE")
    DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
    DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
    DEGenes <- foldChange[names(DEGenes)]
  } else {
    ######## TRADITIONAL WAY TO CHOOSE DE GENES ########
    #print("Classic")
    adj.pvalues <- p.adjust(pvalues, method = 'fdr')
    DEGenes <- foldChange[which(abs(foldChange[adj.pvalues<PVAL]) > 1.5)]
  }
  
  print(length(DEGenes))
}

