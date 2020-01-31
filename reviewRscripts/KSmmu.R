# library(limma)
# y = getGeneKEGGLinks(species.KEGG = "mmu", convert = FALSE)
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/HighEdgeS/inst/extdata/kpgNames_mmu_UpdatedPathways.RData")
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/data/kpg_mmu_UpdatedPathways.RData")
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/data/KEGGgraph_mmu_UpdatedPathways.RData")
library("pROC")

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
j = 1
thresholdFC = 0.5
i = 6

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


temp=data.frame(row.names=names(kpn), Name=kpn,pvalue=rep(NA,length(kpn)))
for (pw in 1:length(kpn)) {
  DEhit <- DEGenes[intersect(keggnodes[[pw]],names(DEGenes))]
  DEmiss <- DEGenes[!DEGenes%in%DEhit]
  
  
  if (length(DEhit) == 0) { pval = NA;}
  else {
    ks <- ks.test(DEhit,DEmiss)
    pval = ks$p.value
  }
  
  temp[names(kpn)[pw],"pvalue"] <- pval
}
temp = temp[order(temp[,2]),]
temp = na.omit(temp)
# myData[[i]] <- temp

mylist <- list()
for (i in 1:length(keggnodes)) {
  mylist[[i]] <- length(intersect(keggnodes[[i]], targetID))
  mylist[[i]]
}
mylist = as.vector(mylist)
names(mylist) = names(keggnodes)
posPredict = mylist[which(mylist == 1)]
negPredict = mylist[!mylist %in% posPredict]

# x <- myData[[1]]
backup <- temp

label <- rep(NA, nrow(temp))
names(label) <- rownames(temp)
label[rownames(temp)%in%names(posPredict)] = 1
label[rownames(temp)%in%names(negPredict)] = 0

temp = cbind.data.frame(temp,label)
colnames(temp) = c("Pathway", "pvalue", "label")
temp = na.omit(temp)
#temp

labels = as.factor(temp[,"label"])


########## IF there is only 1 level in label  #############################
if (nlevels(labels) == 1) {
  temp1 <- cbind.data.frame(Pathway = x$kpn, pvalue = rep(1, length(names(keggnodes))))
  temp1[rownames(temp),2] <- temp[rownames(temp), 2]
  temp1 <- temp1[order(temp1[,2]),]
  
  temp <- temp1
  
  label <- rep(NA, nrow(temp))
  names(label) <- rownames(temp)
  label[rownames(temp)%in%names(posPredict)] = 1
  label[rownames(temp)%in%names(negPredict)] = 0
  
  temp = cbind.data.frame(temp,label)
  colnames(temp) = c("Pathway", "pvalue", "label")
  temp = na.omit(temp)
  labels = as.factor(temp[,"label"])
}

############################################################################
rocobj <- roc(labels, as.numeric(temp[,"pvalue"]),percent = T)  
# pdf(paste0("ROC_KS_", dataset,"_", targetSymbol, "_",methodDE[j],".pdf"), width = 5, height = 5)
# plot(rocobj,legacy.axes=TRUE, main= "ROC curve of Kolmogorov-Smirnov", percent=TRUE, col="red", print.auc=T )
# dev.off()
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


save(x, temp, rocobj, TP, FP, FN, TN, sensitivity, accuracy, specificity,  
     file = paste0("ks",dataset, "_",targetSymbol, "_", methodDE[j], ".RData"))






############ SAMER's DATASETS ###############
load("/Users/GaMinh/Downloads/kegg_mmu.RData")
kpg = x$kpg
kpn = x$kpn
getnodes=function(x){
  return (x@nodes);}
keggnodes=lapply(kpg,getnodes)


filePath = c("/Users/GaMinh/Downloads/altEntrez_GSE22873_Myd88_mmu_SAMER.RData", 
             "Code and Data/data/altEntrez_GSE6030_Neurod1_mmu_SAMER.RData",
             "Code and Data/data/altEntrez_GSE29048_Pdx1_mmu_SAMER.RData")
EsetFile = c("Code and Data/data/DatasetsGAGE/GSE22873_NORM/GSE22873_altCdf.RData",
             "Code and Data/data/DatasetsGAGE/GSE6030_NORM/GSE6030_altCdf.RData",
             "Code and Data/data/DatasetsGAGE/GSE29048_NORM/GSE29048_altCdf.RData")

datasets = c("GSE22873_Myd88","GSE6030NEUROD1", "GSE29048_Pdx1")
targetGene = c("mmu:17874", "mmu:18012", "mmu:18609") 
targetSymbols = c("MYD88", "NEUROD1", "PDX1")

methodDE = c("Classic", "400DE")
j = 1
i = 3
#mmu:17874 = Myd88
#mmu:18012 = NEUROD1
#mmu:18609 = PDX1

load(filePath[i])
load(EsetFile[i])
dataset = datasets[i]
targetID = targetGene[i]
targetSymbol = targetSymbols[i]


# source("https://bioconductor.org/biocLite.R")
# biocLite("BLMA")
# library("BLMA")
# library("ROntoTools")
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

res=data.frame(row.names=names(kpn), Name=kpn,pvalue=rep(NA,length(kpn)))
for (pw in 1:length(kpn)) {
  #print(kpn[pw])
  #pwName <- names(kpn)[pw]
  #DEhit <- DEGenes[intersect(get(paste0("keggnodes$`",pwName,"`")),names(DEGenes))]
  DEhit <- DEGenes[intersect(keggnodes[[pw]],names(DEGenes))]
  DEmiss <- DEGenes[!DEGenes%in%DEhit]
  
  
  if (length(DEhit) == 0) { pval = NA;}
  else {
    ks <- ks.test(DEhit,DEmiss)
    pval = ks$p.value
  }
  
  res[names(kpn)[pw],"pvalue"] <- pval
}
res = res[order(res[,2]),]
temp <- na.omit(res)
# myData[[i]] <- temp

mylist <- list()
for (i in 1:length(keggnodes)) {
  mylist[[i]] <- length(intersect(keggnodes[[i]], targetID))
  mylist[[i]]
}
mylist = as.vector(mylist)
names(mylist) = names(keggnodes)
posPredict = mylist[which(mylist == 1)]
negPredict = mylist[!mylist %in% posPredict]

# x <- myData[[1]]
backup <- temp

label <- rep(NA, nrow(temp))
names(label) <- rownames(temp)
label[rownames(temp)%in%names(posPredict)] = 1
label[rownames(temp)%in%names(negPredict)] = 0

temp = cbind.data.frame(temp,label)
colnames(temp) = c("Pathway", "pvalue", "label")
temp = na.omit(temp)
temp

labels = factor(temp[,"label"])
########## IF there is only 1 level in label  #############################
if (nlevels(labels) ==1 ) {
  temp1 <- cbind.data.frame(Pathway = x$kpn, pvalue = rep(1, length(names(keggnodes))))
  temp1[rownames(temp),2] <- temp[rownames(temp), 2]
  temp1 <- temp1[order(temp1[,2]),]
  
  temp <- temp1
  
  label <- rep(NA, nrow(temp))
  names(label) <- rownames(temp)
  label[rownames(temp)%in%names(posPredict)] = 1
  label[rownames(temp)%in%names(negPredict)] = 0
  
  temp = cbind.data.frame(temp,label)
  colnames(temp) = c("Pathway", "pvalue", "label")
  temp = na.omit(temp)
  labels = factor(temp[,"label"])
} 
print(nlevels(labels))
############################################################################


rocobj <- roc(labels, as.numeric(temp[,"pvalue"]),percent = T)  
pdf(paste0("ROC_KS_", dataset,"_", targetSymbol, "_",methodDE[j],".pdf"), width = 5, height = 5)
plot(rocobj,legacy.axes=TRUE, main= "ROC curve of Kolmogorov-Smirnov", percent=TRUE, col="red", print.auc=T )
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


save(x, temp, rocobj, TP, FP, FN, TN, sensitivity, accuracy, specificity,  
     file = paste0("ks",dataset, "_",targetSymbol, "_", methodDE[j], ".RData"))
