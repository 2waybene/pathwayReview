# install.packages(paste0(baseDir,"BLMA_0.99.4.tar.gz"), repos = NULL,
                 # type="source")
# source("https://bioconductor.org/biocLite.R")
# biocLite("BLMA")
# biocLite("mouse4302.db")


library("pROC")

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

i = 2

dataset = datasets[i]
targetID = targetGenes[i]
targetSymbol = targetGenesSymbol[i]

filePath = paste0(path, "dataset/" , dataset,"_", targetSymbol, "/", dataset, "_", targetSymbol, ".RData")
load(filePath)

data <- get(paste("gene_",dataset,sep=""))
group <- get(paste("group_",dataset,sep=""))
rngroup <- group$Sample
group <- group$Group
names(group) <- rngroup

dataList = list(data)
names(dataList) = dataset
groupList = list(group)
names(groupList) = dataset

GSAComb <- bilevelAnalysisGeneset(gslist, gs.names, splitSize =
                                    ncol(dataList[[1]]), dataList, groupList, enrichment = "GSA", 
                                  nperms = 200, random.seed = 1)

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
temp <- GSAComb
pwID <- as.vector(rownames(temp))
pwName <- as.vector(temp[,"Name"])
pvalue <- as.vector(temp[,"pBLMA"])
temp = cbind.data.frame(pwName,pvalue)
rownames(temp) <- pwID

label <- rep(NA, nrow(temp))
names(label) <- rownames(temp)
label[names(posPredict)] = 1
label[names(negPredict)] = 0

temp = cbind.data.frame(temp,label)
colnames(temp) = c("Pathway", "pvalue", "label")
temp = na.omit(temp)

gsatest = temp

labels = factor(temp[,"label"])
rocobj <- roc(labels, as.numeric(temp[,"pvalue"]),percent = T)  
pdf(paste0("ROC_GSA_", dataset,".pdf"), width = 5, height = 5)
plot(rocobj,legacy.axes=TRUE, main= "ROC curve of GSA", percent=TRUE, col="red", print.auc=T )
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
options(digits=2)
# Calculate sensitivity
sensitivity = TP/(TP+FN)
sensitivity
accuracy = (TN+TP)/(TP+FP+FN+TN)
accuracy
specificity = TN/(TN+FP)
specificity
options(digits=3)
auc(roc(labels, as.numeric(temp[,"pvalue"])))

save(GSAComb, temp, rocobj, TP, FP, FN, TN, sensitivity, accuracy, specificity,  file = paste0("gsa",dataset, "_",targetSymbol, ".RData"))









############ SAMER's DATASETS ###############

library("BLMA")
library("mouse4302.db")
library("pROC")
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/pathwayMouseMapping.RData")
load("/Users/GaMinh/Downloads/kegg_mmu.RData")
kpg = x$kpg
kpn = x$kpn
getnodes=function(x){
  return (x@nodes);}
keggnodes=lapply(kpg,getnodes)


# x <- loadKEGGPathways()
gslist <- lapply(x$kpg,FUN=function(y){return (nodes(y));})
gs.names <- x$kpn[names(gslist)]

path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/mmu/"
setwd(path);
load("/Users/GaMinh/Downloads/kegg_mmu.RData")
i = 1
filePath = c("/Users/GaMinh/Downloads/altEntrez_GSE22873_Myd88_mmu_SAMER.RData", 
             "Code and Data/data/altEntrez_GSE6030_Neurod1_mmu_SAMER.RData",
             "Code and Data/data/altEntrez_GSE29048_Pdx1_mmu_SAMER.RData")
EsetFile = c("Code and Data/data/DatasetsGAGE/GSE22873_NORM/GSE22873_altCdf.RData",
             "Code and Data/data/DatasetsGAGE/GSE6030_NORM/GSE6030_altCdf.RData",
             "Code and Data/data/DatasetsGAGE/GSE29048_NORM/GSE29048_altCdf.RData")

datasets = c("GSE22873_Myd88","GSE6030NEUROD1", "GSE29048_Pdx1")
targetGene = c("mmu:17874", "mmu:18012", "mmu:18609") 
#mmu:17874 = Myd88
#mmu:18012 = NEUROD1
#mmu:18609 = PDX1

load(filePath[i])
load(EsetFile[i])
dataset = datasets[i]
target = targetGene[i]

# create expsMat
samples = c(DataObject$gpControl, DataObject$gpCase)
group = c(rep("c",length(DataObject$gpControl)), rep("d",length(DataObject$gpCase)))
names(group) = samples
mydat = exprs(eset)
mydat = mydat[,samples]

#data_GSE22873 = mydat
#group_GSE22873 = group

#save(data_GSE22873, group_GSE22873, file = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_22/GSE22873.RData")

rownames(mydat) = paste0("mmu:",rownames(mydat))
dataList = list(mydat)
names(dataList) = dataset
groupList = list(group)
names(groupList) = dataset

# perform bi-level meta-analysis in conjunction with ORA
GSAComb <- bilevelAnalysisGeneset(gslist, gs.names, splitSize =
                                    ncol(dataList[[1]]), dataList, groupList, enrichment = "GSA", 
                                  nperms = 200, random.seed = 1)

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
temp <- GSAComb
pwID <- as.vector(rownames(temp))
pwName <- as.vector(temp[,"Name"])
pvalue <- as.vector(temp[,"pBLMA"])
temp = cbind(pwName,pvalue)
rownames(temp) <- pwID

label <- rep(NA, nrow(temp))
names(label) <- rownames(temp)
label[names(posPredict)] = 1
label[names(negPredict)] = 0

temp = cbind.data.frame(temp,label)
colnames(temp) = c("Pathway", "pvalue", "label")
temp = na.omit(temp)

gsatest = temp

labels = factor(temp[,"label"])
rocobj <- roc(labels, as.numeric(temp[,"pvalue"]),percent = T)  
#plot(rocobj,legacy.axes=TRUE, main=plotName, percent=TRUE, col="black", print.auc=T)
plot(rocobj,legacy.axes=TRUE, main= "ROC curve of GSA", percent=TRUE, col="red", print.auc=T )

# Confusion matrix:
TP = nrow(temp[as.numeric(temp[,"pvalue"]) <= 0.05 & temp[,"label"] == 1, ])
FP = nrow(temp[as.numeric(temp[,"pvalue"]) <= 0.05 & temp[,"label"] == 0, ])
FN = nrow(temp[as.numeric(temp[,"pvalue"]) > 0.05 & temp[,"label"] == 1, ])
TN = nrow(temp[as.numeric(temp[,"pvalue"]) > 0.05 & temp[,"label"] == 0, ])

# Calculate sensitivity
sensitivity = TP/(TP+FN)
sensitivity
accuracy = (TN+TP)/(TP+FP+FN+TN)
accuracy
specificity = TN/(TN+FP)
specificity
auc(roc(labels, as.numeric(temp[,"pvalue"])))

#mappingPathway <- cbind(names(kpn), kpn)

save(GSAComb, temp, rocobj, TP, FP, FN, TN, sensitivity, accuracy, specificity,  file = paste0("gsa",dataset, ".RData"))
