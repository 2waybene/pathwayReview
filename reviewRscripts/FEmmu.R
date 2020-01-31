# library(limma)
# y = getGeneKEGGLinks(species.KEGG = "mmu", convert = FALSE)
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/HighEdgeS/inst/extdata/kpgNames_mmu_UpdatedPathways.RData")
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/data/kpg_mmu_UpdatedPathways.RData")
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/data/KEGGgraph_mmu_UpdatedPathways.RData")
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

load("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/mmu/kegg_mmu.RData")
kpg = x$kpg
kpn = x$kpn

getnodes=function(x){
  return (x@nodes);}
keggnodes=lapply(kpg,getnodes)
iterno=2000
PVAL <- 0.05
maxDE=400

path="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/mmu/"
setwd(path);

datasets = c("GSE70302", "GSE70302", "GSE58120", "GSE46211", "GSE49166", "GSE50933", "GSE62999", "GSE57917")
targetGenes = c("mmu:16175", "mmu:16176", "mmu:16428", "mmu:21813", "mmu:20893", "mmu:15903", "mmu:240672", "mmu:15379") 
targetGenesSymbol = c("IL1a", "IL1b", "IL2", "TGFBR2", "BHLHE40", "ID3", "DUSP5", "ONECUR1")

methodDE = c("Classic", "400DE")
j = 2
thresholdFC = 0.5
i = 5 #1, 2, 4, 5

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

length(DEGenes)

####### Saving DE Genes for Sorin #############
# pvalueSymbol <- adj.pvalues
# names(pvalueSymbol) <- vectorMapping[gsub("mmu:","",names(pvalueSymbol))]
# tableIPG <- cbind.data.frame(names(pvalueSymbol), pvalueSymbol, foldChange)
# colnames(tableIPG) <- c("Gene", "Pvalue", "FoldChange")
# write.table(tableIPG, file = paste0(dataset,targetSymbol, "tableIPG.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


###### Done saving DE genes for Sorin ########

tempFE=data.frame(row.names=names(kpn), Name=kpn,pvalue=rep(NA,length(kpn)))
for (pw in 1:length(kpn)) {
  #print(kpn[pw])
  DEinPW = length(intersect(keggnodes[[pw]],names(DEGenes)))
  #print(DEinPW)
  DEnotInPW = length(DEGenes) - DEinPW
  notDEinPW = length(keggnodes[[pw]]) - DEinPW
  notDEnotinPW = length(pvalues) - length(DEGenes) - length(keggnodes[[pw]]) + DEinPW
  table <- matrix(c(DEinPW, DEnotInPW, notDEinPW, notDEnotinPW), nrow = 2, dimnames = list(Genes = c("in PW", "not in PW"),Pathway = c("DE", "not DE")))
  fisher <- fisher.test(table,alternative = "greater")
  #fisher$p.value
  tempFE[pw,"PathwaySize"] = length(keggnodes[[pw]])
  tempFE[pw, "NDE"] = DEinPW
  tempFE[pw,"pvalue"] <- fisher$p.value
}
tempFE = tempFE[order(tempFE[,2]),]
# myData[[i]] <- tempFE

mylist <- lapply(keggnodes, function(x) { length(intersect(x, targetID))} )
mylist = as.vector(mylist)
names(mylist) = names(keggnodes)
posPredict = mylist[which(mylist == 1)]
negPredict = mylist[!mylist %in% posPredict]

# x <- myData[[1]]
# x <- tempFE
# pwID <- rownames(tempFE)
# pwName <- as.vector(tempFE$Name)
# tempFE = cbind.data.frame(pwName,tempFE$pvalue)
# rownames(tempFE) <- pwID

label <- rep(NA, nrow(tempFE))
names(label) <- rownames(tempFE)
label[names(posPredict)] = 1
label[names(negPredict)] = 0

tempFE = cbind.data.frame(tempFE,label)
colnames(tempFE) = c("Name", "pFisher", "PathwaySize", "NDE", "TargetPathway")
tempFE <- tempFE[,c("Name",  "TargetPathway", "PathwaySize", "NDE", "pFisher")]
write.csv(tempFE, file = paste0(datasets[i], "_", targetGenesSymbol[i], ".csv"))

tempFE = na.omit(tempFE)
# tempFE = as.data.frame(tempFE)
fetest <- tempFE

labels = factor(tempFE[,"label"])
rocobj <- roc(labels, as.numeric(tempFE[,"pvalue"]),percent = T)  
#plot(rocobj,legacy.axes=TRUE, main=plotName, percent=TRUE, col="black", print.auc=T)
#pdf(paste0("ROC_FE_", dataset,"_", targetSymbol, "_",methodDE[j],".pdf"), width = 5, height = 5)
#plot(rocobj,legacy.axes=TRUE, main= "ROC curve of Fisher's Exact Test", percent=TRUE, col="red", print.auc=T )
#dev.off()
# Confusion matrix:
TP = nrow(tempFE[as.numeric(tempFE[,"pvalue"]) <= 0.05 & tempFE[,"label"] == 1, ])
FP = nrow(tempFE[as.numeric(tempFE[,"pvalue"]) <= 0.05 & tempFE[,"label"] == 0, ])
FN = nrow(tempFE[as.numeric(tempFE[,"pvalue"]) > 0.05 & tempFE[,"label"] == 1, ])
TN = nrow(tempFE[as.numeric(tempFE[,"pvalue"]) > 0.05 & tempFE[,"label"] == 0, ])
TP
FP
FN
TN
# Calculate sensitivity
options(digits=2)
sensitivity = TP/(TP+FN)
sensitivity
accuracy = (TN+TP)/(TP+FP+FN+TN)
accuracy
specificity = TN/(TN+FP)
specificity
options(digits=3)
auc(roc(labels, as.numeric(tempFE[,"pvalue"])))

save(x, tempFE, rocobj, TP, FP, FN, TN, sensitivity, accuracy, specificity,  
     file = paste0("fe",dataset, "_",targetSymbol,"_",methodDE[j], ".RData"))




############### SAMER's DATA SETS #######################

load("/Users/GaMinh/Downloads/kegg_mmu.RData")
kpg = x$kpg
kpn = x$kpn

getnodes=function(x){
  return (x@nodes);}
keggnodes=lapply(kpg,getnodes)

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
thresholdFC = 1.3
i = 1
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
  print("400DE")
  DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
  DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
  DEGenes <- foldChange[names(DEGenes)]
} else {
  ######## TRADITIONAL WAY TO CHOOSE DE GENES ########
  print("Classic")
  adj.pvalues <- p.adjust(pvalues, method = 'fdr')
  DEGenes <- foldChange[which(abs(foldChange[adj.pvalues<PVAL]) > 1.5)]
}

length(DEGenes)

tempFE=data.frame(row.names=names(kpn), Name=kpn,pvalue=rep(NA,length(kpn)))
for (pw in 1:length(kpn)) {
  #print(kpn[pw])
  DEinPW = length(intersect(keggnodes[[pw]],names(DEGenes)))
  #print(DEinPW)
  DEnotInPW = length(DEGenes) - DEinPW
  notDEinPW = length(keggnodes[[pw]]) - DEinPW
  notDEnotinPW = length(pvalues) - length(DEGenes) - length(keggnodes[[pw]]) + DEinPW
  table <- matrix(c(DEinPW, DEnotInPW, notDEinPW, notDEnotinPW), nrow = 2, dimnames = list(Genes = c("in PW", "not in PW"),Pathway = c("DE", "not DE")))
  fisher <- fisher.test(table)
  #fisher$p.value
  tempFE[pw,"pvalue"] <- fisher$p.value
}
tempFE = tempFE[order(tempFE[,2]),]
# myData[[i]] <- tempFE
mylist <- list()
for (i in 1:length(keggnodes)) {
  mylist[[i]] <- length(intersect(keggnodes[[i]], target))
  mylist[[i]]
}
mylist = as.vector(mylist)
names(mylist) = names(keggnodes)
posPredict = mylist[which(mylist == 1)]
negPredict = mylist[!mylist %in% posPredict]

# x <- myData[[1]]
x <- tempFE
pwID <- rownames(tempFE)
pwName <- as.vector(tempFE$Name)
tempFE = cbind(pwName,tempFE$pvalue)
rownames(tempFE) <- pwID

label <- rep(NA, nrow(tempFE))
names(label) <- rownames(tempFE)
label[names(posPredict)] = 1
label[names(negPredict)] = 0

tempFE = cbind.data.frame(tempFE,label)
colnames(tempFE) = c("Pathway", "pvalue", "label")
tempFE = na.omit(tempFE)
# tempFE = as.data.frame(tempFE)
fetest = tempFE

labels = factor(tempFE[,"label"])
rocobj <- roc(labels, as.numeric(tempFE[,"pvalue"]),percent = T)  
#plot(rocobj,legacy.axes=TRUE, main=plotName, percent=TRUE, col="black", print.auc=T)
pdf(paste0("ROC_FE_", dataset,"_", targetSymbol,"_",methodDE[j],".pdf"), width = 5, height = 5)
plot(rocobj,legacy.axes=TRUE, main= "ROC curve of Fisher's Exact Test", percent=TRUE, col="red", print.auc=T )
dev.off()
# Confusion matrix:
TP = nrow(tempFE[as.numeric(tempFE[,"pvalue"]) <= 0.05 & tempFE[,"label"] == 1, ])
FP = nrow(tempFE[as.numeric(tempFE[,"pvalue"]) <= 0.05 & tempFE[,"label"] == 0, ])
FN = nrow(tempFE[as.numeric(tempFE[,"pvalue"]) > 0.05 & tempFE[,"label"] == 1, ])
TN = nrow(tempFE[as.numeric(tempFE[,"pvalue"]) > 0.05 & tempFE[,"label"] == 0, ])
TP
FP
FN
TN
# Calculate sensitivity
options(digits=2)
sensitivity = TP/(TP+FN)
sensitivity
accuracy = (TN+TP)/(TP+FP+FN+TN)
accuracy
specificity = TN/(TN+FP)
specificity
options(digits=3)
auc(roc(labels, as.numeric(tempFE[,"pvalue"])))

save(x, tempFE, rocobj, TP, FP, FN, TN, sensitivity, accuracy, specificity,  file = paste0("fe",dataset, "_",targetSymbol,"_",methodDE[j], ".RData"))
