##==========================================================
##  Minhs_preprocess_JYL.R
##  Apparently, the file are inheritated rom Minhs
##  All needed libraries are installed properly
##==========================================================

library("hgu133plus2.db")
library("hgu95a.db")
library("hgu133a.db")
library("hgfocus.db")
library("hugene10sttranscriptcluster.db")
library(affyPLM)
library(affy)
library(simpleaffy)
library("GEOquery")
library("ROntoTools")
library("hgu133a2.db")



path="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/"
setwd(path)

load("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/KEGG65.150Pathways.RData")
getnodes=function(x){ return (x@nodes);}
kpg <- setEdgeWeights(kpg)
kpg <- setNodeWeights(kpg, defaultWeight = 1)
keggnodes=lapply(kpg,getnodes)
allGenes=unique(unlist(keggnodes))

design = "Not Paired"
#filename="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/40_new_datasets/datasetslist.txt"
#dataSets = as.character(t(read.table(filename,header=F,sep="\t",stringsAsFactors=F)))
dataset <- "GSE24739"
subgroup <- "G1"
print(dataset)

baseDir <- "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/"
# dataDir <- paste0(baseDir,"dataset/40_new_datasets/",dataset,"/")
dataDir <- paste0(baseDir,"dataset/",dataset, subgroup, "/")
repoDir <- paste0(baseDir,"repository/",dataset)
if (!dir.exists(paste0(baseDir,"repository/",dataset)))
   dir.create(paste0(baseDir,"repository/",dataset), showWarnings = TRUE, recursive = FALSE, mode = "0777")
if (!dir.exists(dataDir))
  dir.create(dataDir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

#filename=paste0(dataDir,dataset,"_annotation.txt")
#annotation=as.character(t(read.table(filename,header=F,sep="\t",stringsAsFactors=F)))
annotation = "hgu133plus2.db"
# gset <- getGEO(filename = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/repository/GSE15824_series_matrix.txt.gz")
gset <- getGEO(dataset, GSEMatrix =TRUE)

if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1 #check the platform
gset <- gset[[idx]]

class <- pData(gset)
head(class, 200)
# temp = class[,c("risk group:ch1", "tissue:ch1")]
# class <- class[class$`race:ch1` == "Caucasian",]
# class <- class["title"]
class <- pData(gset)["title"]
names <- rownames(class)
class <- as.vector(class[,1])
names(class) <- names

xx <- names(class)
head(class, 100)
# class[which(temp$characteristics_ch1 == "cancer type: NSCLC" & temp$characteristics_ch1.1 == "normal: yes")] = "c"
# class[which(temp$characteristics_ch1 == "cancer type: NSCLC" & temp$characteristics_ch1.1 == "normal: no")] = "d"
# class = class[-which(temp$characteristics_ch1.1 == "tissue: medulloblastoma")]
class[grep("G1_Normal",class)]= "c"
class[grep("G1_CML",class)]= "d"
class
class = class[which(class == "c" | class =="d")]
# class[!(class %in% "c")]= "d"
class = as.factor(class)
length(class)
length(xx)

notIncluded <- which(!(xx %in% names(class)))
fileNames <- pData(gset)[-notIncluded,'supplementary_file']
fileNames <- gsub(".*[/]","",fileNames)
rawData <- ReadAffy(verbose=TRUE, celfile.path=repoDir, filenames = fileNames)

# fileNames <- pData(gset)['supplementary_file']
# fileNames$supplementary_file <- gsub(".*[/]","",fileNames$supplementary_file)
# rawData <- ReadAffy(verbose=TRUE, celfile.path=repoDir, filenames = fileNames$supplementary_file)

eset = threestep(rawData,background.method="RMA.2",normalize.method="quantile",summary.method="median.polish")
mydat.d = exprs(eset)
# colnames(mydat.d) <- gsub(".CEL.gz","",colnames(mydat.d))
# mydat.d <- mydat.d[,names(class)]
colnames(mydat.d) = names(class)
temp = data.frame(Sample = names(class), Group = as.vector(class))
rownames(temp) = names(class)
mydat.d.groups <- temp
head(mydat.d)

dynamicVariableData <- paste("data_",dataset, subgroup, sep="")
dynamicVariableGroup <- paste("group_",dataset, subgroup,  sep="")
dynamicVariableAnno <- paste("annotation_",dataset, subgroup, sep="")
dynamicVariableDesign <- paste("design_",dataset, subgroup, sep="")

assign(dynamicVariableData, mydat.d) # Normalized Log
assign(dynamicVariableGroup, mydat.d.groups)
assign(dynamicVariableAnno, annotation)
assign(dynamicVariableDesign, design)

data <- get(paste("data_",dataset,subgroup, sep=""))
group <- get(paste("group_",dataset,subgroup, sep=""))

ano=paste(substring(annotation, 1, nchar(annotation)-3),"ENTREZID",sep="")
require(annotation, character.only = TRUE)
geneAno=as.data.frame(get(ano))
geneAnno=geneAno[,2]
names(geneAnno)=geneAno[,1]
allProbes=unique(geneAno$probe_id)

data <- data[rownames(data)%in%allProbes,]
rownames(data) <- as.character(geneAnno[rownames(data)])
data <- aggregate(data,by=list(rownames(data)),FUN=median)
data <- data[paste("hsa:",data$Group.1,sep="")%in%allGenes,]
rownames(data)=data$Group.1
data <- data[,!colnames(data)%in%c("Group.1")]
rownames(data) = paste0("hsa:",rownames(data))
# colnames(data) = names(class)
head(data)
dynamicVariableGene <- paste("gene_",dataset, subgroup, sep="")
assign(dynamicVariableGene, data) # Normalized Log
range(data)
mylist=c(paste("annotation_",dataset, subgroup, sep=""), 
         paste("design_",dataset, subgroup, sep=""), 
         paste("data_",dataset, subgroup, sep=""), 
         paste("group_",dataset, subgroup, sep=""),
         paste("gene_",dataset, subgroup, sep=""))
save(list=mylist, file=paste0(dataDir,dataset, subgroup, ".RData"))

# setwd("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/40_new_datasets/");
# path="/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/dataset/40_new_datasets/"
# filename="datasetslist.txt"
# dataSets = as.character(t(read.table(filename,header=F,sep="\t",stringsAsFactors=F)))
# for (i in 1:length(dataSets)) {
#   dataset = dataSets[i]
#   print(dataset)
#   dataset = "GSE79962"
#   load(paste0(dataset,"/",dataset, ".RData"))
#   data <- get(paste("data_",dataset,sep=""))
#   print(ncol(data))
#   group <- get(paste("group_",dataset,sep=""))
#   annotation <- get(paste0("annotation_", dataset))
# 
#   ano=paste(substring(annotation, 1, nchar(annotation)-3),"ENTREZID",sep="")
#   require(annotation, character.only = TRUE)
#   geneAno=as.data.frame(get(ano))
#   geneAnno=geneAno[,2]
#   names(geneAnno)=geneAno[,1]
#   allProbes=unique(geneAno$probe_id)
#   
#   data <- data[rownames(data)%in%allProbes,]
#   rownames(data) <- as.character(geneAnno[rownames(data)])
#   data <- aggregate(data,by=list(rownames(data)),FUN=median)
#   data <- data[paste("hsa:",data$Group.1,sep="")%in%allGenes,]
#   rownames(data)=data$Group.1
#   data <- data[,!colnames(data)%in%c("Group.1")]
#   rownames(data) = paste0("hsa:",rownames(data))
#   
#   head(data)
#   dynamicVariableGene <- paste("gene_",dataset, sep="")
#   assign(dynamicVariableGene, data) # Normalized Log
#   range(data)
#   mylist=c(paste("annotation_",dataset, sep=""), 
#            paste("design_",dataset, sep=""), 
#            paste("data_",dataset, sep=""), 
#            paste("group_",dataset, sep=""),
#            paste("gene_",dataset, sep=""))
#   dataDir <- paste0(baseDir,"dataset/40_new_datasets/",dataset,"/")
#   save(list=mylist, file=paste0(dataDir,dataset,".RData"))
# }
