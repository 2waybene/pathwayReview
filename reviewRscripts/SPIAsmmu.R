# library(limma)
# y = getGeneKEGGLinks(species.KEGG = "mmu", convert = FALSE)
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/HighEdgeS/inst/extdata/kpgNames_mmu_UpdatedPathways.RData")
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/data/kpg_mmu_UpdatedPathways.RData")
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/data/KEGGgraph_mmu_UpdatedPathways.RData")
library("pROC")
library(SPIA)
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
names(DEGenes) <- gsub("mmu:","",names(DEGenes))
names(foldChange) <- gsub("mmu:", "", names(foldChange))

res=spia(de=DEGenes, all=names(foldChange), organism="mmu",beta=NULL,nB=2000,plots=FALSE, verbose=TRUE,combine="fisher")
#peRes=pe(de=DEGenes,graphs=x$kpg,ref=names(foldChange),nboot=iterno)
backup <- res
res <- res[,-12]
rownames(res) <- res[,1]
res[,1] <- NULL
res[,dim(res)[2]+1] <- res[,1]
res[,dim(res)[2]] <- paste0("path:mmu",res[,dim(res)[2]])
res[,1] <- NULL
colnames(res)[dim(res)[2]] = "pathway"
# res

########## IF there are 2 levels in label already ########################
temp <- cbind.data.frame(rownames(res), res$pG)
colnames(temp) <- c("Pathway", "pvalue")
rownames(temp) <- res$pathway
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
temp = na.omit(temp)

spiaRes = temp

########## IF there is only 1 level in label  ########################

mappingPathway <- cbind.data.frame(names(kpn), kpn)
pwname = as.vector(mappingPathway[,1])
names(pwname) = mappingPathway[,2]
rnamesRes = pwname[rownames(res)]
NAs = which(is.na(rnamesRes))
rnamesRes = na.omit(rnamesRes)
res = res[-NAs,]
temp <- cbind.data.frame(rownames(res), res$pG)
rownames(temp) <- rnamesRes
colnames(temp) <- c("pathway", "pvalue")
rnamesTemp = rownames(temp)[complete.cases(rownames(temp))]
temp = temp[rnamesTemp,]

temp1 <- cbind.data.frame(Pathway = x$kpn, pvalue = rep(1, length(names(keggnodes))))
temp1[rownames(temp),2] <- temp[rownames(temp), 2]
temp1 <- temp1[order(temp1[,2]),]

temp <- temp1

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
#colnames(temp) = c("Pathway", "pvalue", "label")
temp = na.omit(temp)

spiaRes = temp

##############################################################

labels = factor(temp[,"label"])
rocobj <- roc(labels, as.numeric(temp[,"pvalue"]),percent = T)  
# pdf(paste0("ROC_SPIA_", dataset,"_", targetSymbol, "_",methodDE[j],".pdf"), width = 5, height = 5)
# plot(rocobj,legacy.axes=TRUE, main= "ROC curve of SPIA", percent=TRUE, col="red", print.auc=T )
# dev.off()

# Confusion matrix:
TP = nrow(temp[as.numeric(temp[,"pvalue"]) < 0.05 & temp[,"label"] == 1, ])
FP = nrow(temp[as.numeric(temp[,"pvalue"]) < 0.05 & temp[,"label"] == 0, ])
FN = nrow(temp[as.numeric(temp[,"pvalue"]) > 0.05 & temp[,"label"] == 1, ])
TN = nrow(temp[as.numeric(temp[,"pvalue"]) > 0.05 & temp[,"label"] == 0, ])
TP
FP
FN
TN
options(digits=2)
# Calculate sensitivity
sensitivity = TP/(TP+FN)
sensitivity
accuracy = (TN+TP)/(TP+FP+FN+TN)
accuracy
specificity = TN/(TN+FP)
specificity
options(digits=3)
auc(roc(labels,as.numeric(temp[,"pvalue"])))
#mappingPathway <- cbind(names(kpn), kpn)

save(temp, spiaRes, rocobj, TP, FP, FN, TN, sensitivity, accuracy, specificity,  
     file = paste0("spia", dataset, "_",targetSymbol,"_",methodDE[1], ".RData"))



############ SAMER's DATASETS ###############
library(SPIA)
library("pROC")
path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/mmu/"
setwd(path);
load("/Users/GaMinh/Downloads/kegg_mmu.RData")

filePath = c("/Users/GaMinh/Downloads/altEntrez_GSE22873_Myd88_mmu_SAMER.RData", 
             "Code and Data/data/altEntrez_GSE6030_Neurod1_mmu_SAMER.RData",
             "Code and Data/data/altEntrez_GSE29048_Pdx1_mmu_SAMER.RData")
datasets = c("GSE22873_Myd88","GSE6030NEUROD1", "GSE29048_Pdx1")
targetGene = c("mmu:17874", "mmu:18012", "mmu:18609") 
targetSymbols = c("MYD88", "NEUROD1", "PDX1")
methodDE = c("Classic", "400DE")

j = 1
i = 3
load(filePath[i])
dataset = datasets[i]
target = targetGene[i]
targetSymbol = targetSymbols[i]

kpg = x$kpg
kpn = x$kpn
getnodes=function(x){
  return (x@nodes);}
keggnodes=lapply(kpg,getnodes)

#library("ROntoTools")
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
names(DEGenes) = gsub("mmu:", "", names(DEGenes))
names(foldChange) = gsub("mmu:", "", names(foldChange))

res=spia(de=DEGenes, all=names(foldChange), organism="mmu",beta=NULL,nB=2000,plots=FALSE, verbose=TRUE,combine="fisher")
#peRes=pe(de=DEGenes,graphs=x$kpg,ref=names(foldChange),nboot=iterno)
res <- res[,-12]
rownames(res) <- res[,1]
res[,1] <- NULL
res[,dim(res)[2]+1] <- res[,1]
res[,dim(res)[2]] <- paste0("path:mmu",res[,dim(res)[2]])
res[,1] <- NULL
colnames(res)[dim(res)[2]] = "pathway"
res

res <- cbind.data.frame(res, rownames(res))
rownames(res) <- res$pathway

temp <- res[c(11,6)]
mylist <- list()
for (i in 1:length(keggnodes)) {
  mylist[[i]] <- length(intersect(keggnodes[[i]], target))
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
temp = na.omit(temp)

# mappingPathway <- cbind(names(kpn), kpn)
# pwname = mappingPathway[,1]
# names(pwname) = mappingPathway[,2]
# rnamesRes = pwname[rownames(res)]
# NAs = which(is.na(rnamesRes))
# rnamesRes = na.omit(rnamesRes)
# res = res[-NAs,]
# temp <- cbind.data.frame(rownames(res), res$pG)
# rownames(temp) <- rnamesRes
# colnames(temp) <- c("pathway", "pvalue")
# rnamesTemp = rownames(temp)[complete.cases(rownames(temp))]
# temp = temp[rnamesTemp,]
###############
# pathwaymapping = cbind(rownames(temp), temp$pathway)
#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/pathwayMouseMapping.RData")
spiaRes = temp

labels = factor(temp[,"label"])

rocobj <- roc(labels, as.numeric(temp[,"pvalue"]),percent = T)  
pdf(paste0("ROC_SPIA_", dataset,"_", targetSymbol, "_",methodDE[j],".pdf"), width = 5, height = 5)
plot(rocobj,legacy.axes=TRUE, main= "ROC curve of SPIA", percent=TRUE, col="red", print.auc=T )
dev.off()
auc(rocobj)

# Confusion matrix:
TP = nrow(temp[as.numeric(temp[,"pvalue"]) < 0.05 & temp[,"label"] == 1, ])
FP = nrow(temp[as.numeric(temp[,"pvalue"]) < 0.05 & temp[,"label"] == 0, ])
FN = nrow(temp[as.numeric(temp[,"pvalue"]) > 0.05 & temp[,"label"] == 1, ])
TN = nrow(temp[as.numeric(temp[,"pvalue"]) > 0.05 & temp[,"label"] == 0, ])

# Calculate sensitivity
sensitivity = TP/(TP+FN)
sensitivity
accuracy = (TN+TP)/(TP+FP+FN+TN)
accuracy
specificity = TN/(TN+FP)
specificity
auc(roc(labels,as.numeric(temp[,"pvalue"])))
#mappingPathway <- cbind(names(kpn), kpn)

save(temp, res, rocobj, TP, FP, FN, TN, sensitivity, accuracy, specificity,  
     file = paste0("spia", dataset, "_",targetSymbol,"_",methodDE[1], ".RData"))


