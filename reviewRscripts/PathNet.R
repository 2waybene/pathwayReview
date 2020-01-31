# source("http://bioconductor.org/biocLite.R")
# biocLite("PathNet")
# biocLite("PathNetData")

library("PathNet")
library("PathNetData")
data(A)
data(pathway)

# load the matrix A for pathnet
# setwd(system.file(dir="extdata", package="PathNetData"))
# A <- as.matrix(read.table(file = "adjacency_data.txt", sep = "\t", header = T))
# pathway <- read.table(file = "pathway_data.txt", sep = "\t", header = T)

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
for (i in 1:length(dataSets)) {
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
                     DirectEvidence_info = pvalueMatrix,
                     Adjacency = A,
                     pathway = pathway,
                     Column_DirectEvidence = 2,
                     n_perm = 2000, threshold = 0.05)
  
  myData[[i]] <- results$enrichment_results
}
names(myData) <- dataSets
#2nd run takes 12000 row of pvalueMatrix
#3nd run takes everything in pvalueMatrix
save(myData,file=paste(path,"alldatasetsPathNet (3nd run).RData",sep=""))

################ number of significant pathway ###########

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsPathNet (3rd run).RData")
signPW = vector(mode = "numeric", length = 0)
for (i in 1:length(myData)) {
  signPW = c(signPW, length(which(myData[[i]][,'p_PathNet']<0.05))/nrow(myData[[i]]))
}

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/PathNet40Datasets.RData")
for (i in 1:length(myData)) {
  signPW = c(signPW, length(which(myData[[i]][,'p_PathNet']<0.05))/nrow(myData[[i]]))
}

mean(signPW)

########## Ranks and Pvalues of target Pathways ########## 
load≈
load("diseaseList.RData")
rankPathNetWithoutNull <- rep(NA, 35)
pValuePathNetWithoutNull <- rep(NA, 35)

for(i in 1:length(Renal_cell_carcinoma)) {
  rankPathNetWithoutNull[Renal_cell_carcinoma[i]] <- rownames(myData[[Renal_cell_carcinoma[i]]][myData[[Renal_cell_carcinoma[i]]][,1]=="Renal cell carcinoma",])
  pValuePathNetWithoutNull[Renal_cell_carcinoma[i]] <- myData[[Renal_cell_carcinoma[i]]][myData[[Renal_cell_carcinoma[i]]][,1]=="Renal cell carcinoma","p_PathNet"]
}

for(i in 1:length(Alzheimers_disease)) {
  rankPathNetWithoutNull[Alzheimers_disease[i]] <- rownames(myData[[Alzheimers_disease[i]]][myData[[Alzheimers_disease[i]]][,1]=="Alzheimers disease",])
  pValuePathNetWithoutNull[Alzheimers_disease[i]] <- myData[[Alzheimers_disease[i]]][myData[[Alzheimers_disease[i]]][,1]=="Alzheimers disease","p_PathNet"]
}

for(i in 1:length(Thyroid_cancer)) {
  rankPathNetWithoutNull[Thyroid_cancer[i]] <- rownames(myData[[Thyroid_cancer[i]]][myData[[Thyroid_cancer[i]]][,1]=="Thyroid cancer",])
  pValuePathNetWithoutNull[Thyroid_cancer[i]] <- myData[[Thyroid_cancer[i]]][myData[[Thyroid_cancer[i]]][,1]=="Thyroid cancer","p_PathNet"]
}

for(i in 1:length(Colorectal_cancer)) {
  rankPathNetWithoutNull[Colorectal_cancer[i]] <- rownames(myData[[Colorectal_cancer[i]]][myData[[Colorectal_cancer[i]]][,1]=="Colorectal cancer",])
  pValuePathNetWithoutNull[Colorectal_cancer[i]] <- myData[[Colorectal_cancer[i]]][myData[[Colorectal_cancer[i]]][,1]=="Colorectal cancer","p_PathNet"]
}

for(i in 1:length(Prostate_cancer)) {
  rankPathNetWithoutNull[Prostate_cancer[i]] <- rownames(myData[[Prostate_cancer[i]]][myData[[Prostate_cancer[i]]][,1]=="Prostate cancer",])
  pValuePathNetWithoutNull[Prostate_cancer[i]] <- myData[[Prostate_cancer[i]]][myData[[Prostate_cancer[i]]][,1]=="Prostate cancer","p_PathNet"]
}

for(i in 1:length(Acute_myeloid_leukemia)) {
  rankPathNetWithoutNull[Acute_myeloid_leukemia[i]] <- rownames(myData[[Acute_myeloid_leukemia[i]]][myData[[Acute_myeloid_leukemia[i]]][,1]=="Acute myeloid leukemia",])
  pValuePathNetWithoutNull[Acute_myeloid_leukemia[i]] <- myData[[Acute_myeloid_leukemia[i]]][myData[[Acute_myeloid_leukemia[i]]][,1]=="Acute myeloid leukemia","p_PathNet"]
}

for(i in 1:length(Pancreatic_cancer)) {
  rankPathNetWithoutNull[Pancreatic_cancer[i]] <- rownames(myData[[Pancreatic_cancer[i]]][myData[[Pancreatic_cancer[i]]][,1]=="Pancreatic cancer",])
  pValuePathNetWithoutNull[Pancreatic_cancer[i]] <- myData[[Pancreatic_cancer[i]]][myData[[Pancreatic_cancer[i]]][,1]=="Pancreatic cancer","p_PathNet"]
}

for(i in 1:length(Non_small_cell_lung_cancer)) {
  rankPathNetWithoutNull[Non_small_cell_lung_cancer[i]] <- rownames(myData[[Non_small_cell_lung_cancer[i]]][myData[[Non_small_cell_lung_cancer[i]]][,1]=="Non-small cell lung cancer",])
  pValuePathNetWithoutNull[Non_small_cell_lung_cancer[i]] <- myData[[Non_small_cell_lung_cancer[i]]][myData[[Non_small_cell_lung_cancer[i]]][,1]=="Non-small cell lung cancer","p_PathNet"]
}

for(i in 1:length(Glioma)) {
  rankPathNetWithoutNull[Glioma[i]] <- rownames(myData[[Glioma[i]]][myData[[Glioma[i]]][,1]=="Glioma",])
  pValuePathNetWithoutNull[Glioma[i]] <- myData[[Glioma[i]]][myData[[Glioma[i]]][,1]=="Glioma","p_PathNet"]
}

for(i in 1:length(Parkinsons_disease)) {
  rankPathNetWithoutNull[Parkinsons_disease[i]] <- rownames(myData[[Parkinsons_disease[i]]][myData[[Parkinsons_disease[i]]][,1]=="Parkinsons disease",])
  pValuePathNetWithoutNull[Parkinsons_disease[i]] <- myData[[Parkinsons_disease[i]]][myData[[Parkinsons_disease[i]]][,1]=="Parkinsons disease","p_PathNet"]
}

for(i in 1:length(Chronic_myeloid_leukemia)) {
  rankPathNetWithoutNull[Chronic_myeloid_leukemia[i]] <- rownames(myData[[Parkinsons_disease[i]]][myData[[Parkinsons_disease[i]]][,1]=="Parkinsons disease",])
  pValuePathNetWithoutNull[Chronic_myeloid_leukemia[i]] <- myData[[Parkinsons_disease[i]]][myData[[Parkinsons_disease[i]]][,1]=="Parkinsons disease","p_PathNet"]
}

for(i in 1:length(Type_II_diabetes_mellitus)) {
  rankPathNetWithoutNull[Type_II_diabetes_mellitus[i]] <- rownames(myData[[Type_II_diabetes_mellitus[i]]][myData[[Type_II_diabetes_mellitus[i]]][,1]=="Type II diabetes mellitus",])
  pValuePathNetWithoutNull[Type_II_diabetes_mellitus[i]] <- myData[[Type_II_diabetes_mellitus[i]]][myData[[Type_II_diabetes_mellitus[i]]][,1]=="Type II diabetes mellitus","p_PathNet"]
}

for(i in 1:length(Dilated_cardiomyopathy)) {
  rankPathNetWithoutNull[Dilated_cardiomyopathy[i]] <- rownames(myData[[Dilated_cardiomyopathy[i]]][myData[[Dilated_cardiomyopathy[i]]][,1]=="Dilated cardiomyopathy",])
  pValuePathNetWithoutNull[Dilated_cardiomyopathy[i]] <- myData[[Dilated_cardiomyopathy[i]]][myData[[Dilated_cardiomyopathy[i]]][,1]=="Dilated cardiomyopathy","p_PathNet"]
}

for(i in 1:length(Huntingtons_disease)) {
  rankPathNetWithoutNull[Huntingtons_disease[i]] <- rownames(myData[[Huntingtons_disease[i]]][myData[[Huntingtons_disease[i]]][,1]=="Huntingtons disease",])
  pValuePathNetWithoutNull[Huntingtons_disease[i]] <- myData[[Huntingtons_disease[i]]][myData[[Huntingtons_disease[i]]][,1]=="Huntingtons disease","p_PathNet"]
}

for(i in 1:length(Endometrial_cancer)) {
  rankPathNetWithoutNull[Endometrial_cancer[i]] <- rownames(myData[[Endometrial_cancer[i]]][myData[[Endometrial_cancer[i]]][,1]=="Endometrial cancer",])
  pValuePathNetWithoutNull[Endometrial_cancer[i]] <- myData[[Endometrial_cancer[i]]][myData[[Endometrial_cancer[i]]][,1]=="Endometrial cancer","p_PathNet"]
}

save(rankPathNetWithoutNull, pValuePathNetWithoutNull, file=paste0(path,"PathNetWithoutNull (3rd run).RData"))

boxplot(as.vector(as.numeric(rankPathNetWithoutNull)),main ="Ranks of target pathways", ylim=c(0, 130))
boxplot(as.vector(as.numeric(pValuePathNetWithoutNull)),main ="p-values of target pathways", ylim = c(0,1))
library(vioplot)
vioplot(as.vector(as.numeric(rankPathNetWithoutNull)), col = "grey", ylim = c(0,120))
title("Rank of target pathway using PathNet", ylab = "Rank of target pathways")
median(as.vector(as.numeric(rankPathNetWithoutNull)))
vioplot(as.vector(as.numeric(pValuePathNetWithoutNull)), col = "grey", ylim = c(0,1))
title("pValue of target pathway using PathNet", ylab = "Rank of target pathways")
median(as.vector(as.numeric(pValuePathNetWithoutNull)))




