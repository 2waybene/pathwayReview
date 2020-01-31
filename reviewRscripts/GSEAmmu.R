library(clusterProfiler)
library("pROC")
library("ROntoTools")
load("/Users/GaMinh/Downloads/kegg_mmu.RData")
kpg = x$kpg
kpn = x$kpn
getnodes=function(x){
  return (x@nodes);}
keggnodes=lapply(kpg,getnodes)


# x <- loadKEGGPathways()
gslist <- lapply(x$kpg,FUN=function(y){return (nodes(y));})
gs.names <- x$kpn[names(gslist)]

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
j = 1

i = 8

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

names(foldChange) <- gsub("mmu:","", names(foldChange))

foldChange[order(foldChange, decreasing = TRUE)]
res <- gseKEGG(foldChange[order(foldChange, decreasing = TRUE)], organism = 'mmu', keyType = "kegg", 
               exponent = 1, nPerm = 10000, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, 
               pAdjustMethod = "BH", verbose = TRUE, use_internal_data = FALSE, seed = FALSE)

temp = as.data.frame(res)[,c(2,7)]
rownames(temp) = paste0("path:", rownames(temp))

mylist <- list()
for (i in 1:length(keggnodes)) {
  mylist[[i]] <- length(intersect(keggnodes[[i]], targetID))
  mylist[[i]]
}
mylist = as.vector(mylist)
names(mylist) = names(keggnodes)
posPredict = mylist[which(mylist == 1)]
negPredict = mylist[!mylist %in% posPredict]


label <- rep(NA, nrow(temp))
names(label) <- rownames(temp)
label[rownames(temp)%in%names(posPredict)] = 1
label[rownames(temp)%in%names(negPredict)] = 0

temp = cbind.data.frame(temp,label)
colnames(temp) = c("Pathway", "pvalue", "label")
temp <- na.omit(temp)
dim(temp)

labels = factor(temp[,"label"])
rocobj <- roc(labels, as.numeric(temp[,"pvalue"]),percent = T)  
#plot(rocobj,legacy.axes=TRUE, main=plotName, percent=TRUE, col="black", print.auc=T)
pdf(paste0("ROC_GSEA_", dataset,"_", targetSymbol,".pdf"), width = 5, height = 5)
plot(rocobj,legacy.axes=TRUE, main= "ROC curve of GSEA", percent=TRUE, col="red", print.auc=T )
dev.off()
# Confusion matrix:
TP = nrow(temp[as.numeric(temp[,"pvalue"]) <= 0.05 & temp[,"label"] == 1, ])
FP = nrow(temp[as.numeric(temp[,"pvalue"]) <= 0.05 & temp[,"label"] == 0, ])
FN = nrow(temp[as.numeric(temp[,"pvalue"]) > 0.05 & temp[,"label"] == 1, ])
TN = nrow(temp[as.numeric(temp[,"pvalue"]) > 0.05 & temp[,"label"] == 0, ])
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
auc(roc(labels, as.numeric(temp[,"pvalue"])))

save(x, temp, rocobj, TP, FP, FN, TN, sensitivity, accuracy, specificity,  file = paste0("GSEA",dataset,"_",targetSymbol, ".RData"))


############# SAMER's DATA SETS ##################

path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/mmu/"
setwd(path);
load("/Users/GaMinh/Downloads/kegg_mmu.RData")
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

filePath = c("/Users/GaMinh/Downloads/altEntrez_GSE22873_Myd88_mmu_SAMER.RData", 
             "Code and Data/data/altEntrez_GSE6030_Neurod1_mmu_SAMER.RData",
             "Code and Data/data/altEntrez_GSE29048_Pdx1_mmu_SAMER.RData")
datasets = c("GSE22873_Myd88","GSE6030NEUROD1", "GSE29048_Pdx1")
targetGenes = c("mmu:17874", "mmu:18012", "mmu:18609") 
targetGenesSymbol = c("MYD88", "NEUROD1", "PDX1")

i = 3

load(filePath[i])
dataset = datasets[i]
targetID = targetGenes[i]
targetSymbol = targetGenesSymbol[i]

mydata = DataObject$exprTable
rownames(mydata) = paste0("mmu:",rownames(mydata))
foldChange = mydata$logFC
names(foldChange) = gsub("mmu:", "", rownames(mydata))

foldChange[order(foldChange, decreasing = TRUE)]
res <- gseKEGG(foldChange[order(foldChange, decreasing = TRUE)], organism = 'mmu', keyType = "kegg", 
               exponent = 1, nPerm = 10000, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, 
               pAdjustMethod = "BH", verbose = TRUE, use_internal_data = FALSE, seed = FALSE)

temp = as.data.frame(res)[,c(2,7)]
rownames(temp) = paste0("path:", rownames(temp))

mylist <- list()
for (i in 1:length(keggnodes)) {
  mylist[[i]] <- length(intersect(keggnodes[[i]], targetID))
  mylist[[i]]
}
mylist = as.vector(mylist)
names(mylist) = names(keggnodes)
posPredict = mylist[which(mylist == 1)]
negPredict = mylist[!mylist %in% posPredict]


label <- rep(NA, nrow(temp))
names(label) <- rownames(temp)
label[rownames(temp)%in%names(posPredict)] = 1
label[rownames(temp)%in%names(negPredict)] = 0

temp = cbind.data.frame(temp,label)
colnames(temp) = c("Pathway", "pvalue", "label")
temp <- na.omit(temp)
dim(temp)

labels = factor(temp[,"label"])
rocobj <- roc(labels, as.numeric(temp[,"pvalue"]),percent = T)  
#plot(rocobj,legacy.axes=TRUE, main=plotName, percent=TRUE, col="black", print.auc=T)
pdf(paste0("ROC_GSEA_", dataset,"_", targetSymbol,".pdf"), width = 5, height = 5)
plot(rocobj,legacy.axes=TRUE, main= "ROC curve of GSEA", percent=TRUE, col="red", print.auc=T )
dev.off()
# Confusion matrix:
TP = nrow(temp[as.numeric(temp[,"pvalue"]) <= 0.05 & temp[,"label"] == 1, ])
FP = nrow(temp[as.numeric(temp[,"pvalue"]) <= 0.05 & temp[,"label"] == 0, ])
FN = nrow(temp[as.numeric(temp[,"pvalue"]) > 0.05 & temp[,"label"] == 1, ])
TN = nrow(temp[as.numeric(temp[,"pvalue"]) > 0.05 & temp[,"label"] == 0, ])
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
auc(roc(labels, as.numeric(temp[,"pvalue"])))

save(x, temp, rocobj, TP, FP, FN, TN, sensitivity, accuracy, specificity,  file = paste0("GSEA",dataset,"_",targetSymbol, ".RData"))

