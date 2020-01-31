# library(limma)
# y = getGeneKEGGLinks(species.KEGG = "mmu", convert = FALSE)
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/HighEdgeS/inst/extdata/kpgNames_mmu_UpdatedPathways.RData")
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/data/kpg_mmu_UpdatedPathways.RData")
# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/Code and Data/data/KEGGgraph_mmu_UpdatedPathways.RData")
library("pROC")
library("ROntoTools")
load("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/mmu/kegg_mmu.RData")
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

path="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/mmu/"
setwd(path);

datasets = c("GSE70302", "GSE70302", "GSE58120", "GSE46211", "GSE49166", "GSE50933", "GSE62999", "GSE57917")
targetGenes = c("mmu:16175", "mmu:16176", "mmu:16428", "mmu:21813", "mmu:20893", "mmu:15903", "mmu:240672", "mmu:15379") 
targetGenesSymbol = c("IL1a", "IL1b", "IL2", "TGFBR2", "BHLHE40", "ID3", "DUSP5", "ONECUR1")
maxDE=400
methodDE = c("Classic", "400DE")

j = 2
thresholdFC = 0.5
i = 4 #1, 2, 4, 5

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

peRes=pe(x=DEGenes,graphs=x$kpg,ref=names(foldChange),nboot=iterno)

temp <- Summary(peRes)
temp <- cbind.data.frame(temp,x$kpn[rownames(temp)])
colnames(temp) <- c(colnames(Summary(peRes)), "pathway")
# XX <- rep(NA,length(kpn))
# temp1=data.frame(row.names=names(kpn), Name=kpn, pPert=XX, 
#                  pPert.adj=XX, )

#temp <- cbind.data.frame(temp$pathway, temp$pComb)
# rownames(temp) <- rownames(Summary(peRes))
# temp = na.omit(temp)
# head(temp)

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
temp = temp[order(temp$pathway),]
# colnames(temp) = c("Pathway", "pvalue", "label")

write.csv(temp, file = paste0(datasets[i], "_", targetSymbol, "_ROntoTools.csv"))



labels = factor(temp[,"label"]) 

if(nlevels(labels) == 1) {
  ########## IF there is only 1 level in label  #############################
  
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
  
  labels = factor(temp[,"label"]) 
}

rocobj <- roc(labels, as.numeric(temp[,"pvalue"]),percent = T)  
# pdf(paste0("ROC_RONTOTOOLS_", dataset,"_", targetSymbol, "_",methodDE[j],".pdf"), width = 5, height = 5)
# plot(rocobj,legacy.axes=TRUE, main= "ROC curve of R Onto Tools", percent=TRUE, col="red", print.auc=T )
# dev.off()
auc(rocobj)

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

save(peRes, temp, rocobj, TP, FP, FN, TN, sensitivity, accuracy, specificity,  
     file = paste0("rontotools",dataset, "_",targetSymbol,"_",methodDE[j], ".RData"))





############ SAMER's DATASETS ###############

path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/mmu/"
setwd(path);
load("/Users/GaMinh/Downloads/kegg_mmu.RData")
library("pROC")
# pathwaymapping = cbind(rownames(temp), temp$pathway)
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/pathwayMouseMapping.RData")
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
targetSymbols = c("MYD88", "NEUROD1", "PDX1")
targetGene = c("mmu:17874", "mmu:18012", "mmu:18609") 

methodDE = c("Classic", "400DE")
j = 2
i = 1
load(filePath[i])
dataset = datasets[i]
targetID = targetGene[i]
targetSymbol = targetSymbols[i]

mydata = DataObject$exprTable
rownames(mydata) = paste0("mmu:",rownames(mydata))
foldChange = mydata$logFC
names(foldChange) = rownames(mydata)
pvalues = mydata$P.Value
#pvalues = mydata$P.Value
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

peRes=pe(x=DEGenes,graphs=x$kpg,ref=names(foldChange),nboot=iterno)
temp <- Summary(peRes)
temp[,dim(temp)[2]+1] = rownames(temp)
rownames(temp)<-x$kpn[rownames(temp)]
colnames(temp)[dim(temp)[2]] = "pathway"
head(temp)

###############

mylist <- list()
for (i in 1:length(keggnodes)) {
  mylist[[i]] <- length(intersect(keggnodes[[i]], targetID))
  mylist[[i]]
}
mylist = as.vector(mylist)
names(mylist) = names(keggnodes)
posPredict = mylist[which(mylist == 1)]
negPredict = mylist[!mylist %in% posPredict]

#x <- myData[[1]]
x <- temp
pwID <- temp$pathway
pwName <- rownames(temp)
temp = cbind(pwName,temp$pComb)
rownames(temp) <- pwID

label <- rep(NA, nrow(temp))
names(label) <- rownames(temp)
label[rownames(temp)%in%names(posPredict)] = 1
label[rownames(temp)%in%names(negPredict)] = 0

temp = cbind.data.frame(temp,label)
colnames(temp) = c("Pathway", "pvalue", "label")
temp = na.omit(temp)

rontotools = temp

labels = factor(temp[,"label"])
#rocobj <- roc(as.numeric(temp[,"pvalue"]), factor(temp[,"label"]))  
rocobj <- roc(labels, as.numeric(temp[,"pvalue"]),percent = T)  
pdf(paste0("ROC_RONTOTOOLS_", dataset,"_", targetSymbol, "_",methodDE[j],".pdf"), width = 5, height = 5)
plot(rocobj,legacy.axes=TRUE, main= "ROC curve of R Onto Tools", percent=TRUE, col="red", print.auc=T )
dev.off()
auc(rocobj)

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

save(peRes, temp, rocobj, TP, FP, FN, TN, sensitivity, accuracy, specificity,  
     file = paste0("rontotools",dataset, "_",targetSymbol,"_",methodDE[j], ".RData"))



