load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ResultSpia.RData")
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ResultSpiaNormalDEGenes.RData")

rank <- NULL
pValue <- NULL
dataset <- resList$GSE781
disease = "Renal cell carcinoma"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- as.numeric(rownames(dataset[dataset$Name==disease,]))
dataset[dataset$Name==disease, 'pG']
pValue <- dataset[dataset$Name==disease, 'pG']

dataset <- resList$GSE1297
disease = "Alzheimer's disease"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE3467
disease = "Thyroid cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE3585
disease = "Dilated cardiomyopathy"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE3678
disease = "Thyroid cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE4107
disease = "Colorectal cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE5281_EC
disease = "Alzheimer's disease"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE5281_HIP
disease = "Alzheimer's disease"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE5281_VCX
disease = "Alzheimer's disease"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE6956_AA
disease = "Prostate cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE6956_C
disease = "Prostate cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE8671
disease = "Colorectal cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE8762
disease = "Hungtington's disease"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE9348
disease = "Colorectal cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE9476
disease = "Acute myeloid leukemia"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE14762
disease = "Renal cell carcinoma"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE15471
disease = "Pancreatic cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE16515
disease = "Pancreatic cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE18842
disease = "Non-small cell lung cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE19188
disease = "Non-small cell lung cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE19728
disease = "Glioma"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE20153
disease = "Parkinson's disease"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE20291
disease = "Parkinson's disease"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE21354
disease = "Glioma"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE14924_CD4
disease = "Acute myeloid leukemia"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE14924_CD8
disease = "Acute myeloid leukemia"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE16759
disease = "Alzheimer's disease"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE19420
disease = "Type II diabetes mellitus"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE20164
disease = "Parkinson's disease"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE23878
disease = "Colorectal cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE24739_G0
disease = "Chronic myeloid leukemia"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE24739_G1
disease = "Chronic myeloid leukemia"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE32676
disease = "Pancreatic cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE4183
disease = "Colorectal cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

dataset <- resList$GSE7305
disease = "Endometrial cancer"
as.numeric(rownames(dataset[dataset$Name==disease,]))
rank <- cbind(rank,as.numeric(rownames(dataset[dataset$Name==disease,])))
dataset[dataset$Name==disease, 'pG']
pValue <- cbind(pValue,dataset[dataset$Name==disease, 'pG'])

boxplot(as.vector(rank),main ="Ranks of target pathways", ylim=c(0, 120))
boxplot(as.vector(pValue),main ="p-values of target pathways")
hist(rank)
hist(pValue, breaks = 50)
resList$GSE7305
hist(resList$GSE7305)
