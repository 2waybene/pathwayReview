load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ResultSpia.RData")
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ResultSpiaNormalDEGenes.RData")
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/resList_SPIA_400DE.RData")

rank <- NULL
pValue <- NULL
dataset <- resList[[1]]
disease = "Renal cell carcinoma"
match(disease,rownames(dataset))
rank <- match(disease,rownames(dataset))
#rank <- match(disease,rownames(dataset))
dataset[disease, 'pG']
pValue <- dataset[disease, 'pG']

dataset <- resList[[2]]
disease = "Alzheimer's disease"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[3]]
disease = "Thyroid cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[4]]
disease = "Dilated cardiomyopathy"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[5]]
disease = "Thyroid cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[6]]
disease = "Colorectal cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[7]]
disease = "Alzheimer's disease"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[8]]
disease = "Alzheimer's disease"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[9]]
disease = "Alzheimer's disease"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[10]]
disease = "Prostate cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[11]]
disease = "Prostate cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[12]]
disease = "Colorectal cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[13]]
disease = "Hungtington's disease"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[14]]
disease = "Colorectal cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[15]]
disease = "Acute myeloid leukemia"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[16]]
disease = "Renal cell carcinoma"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[17]]
disease = "Pancreatic cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[18]]
disease = "Pancreatic cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[19]]
disease = "Non-small cell lung cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[20]]
disease = "Non-small cell lung cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[21]]
disease = "Glioma"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[22]]
disease = "Parkinson's disease"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[23]]
disease = "Parkinson's disease"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[24]]
disease = "Glioma"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[25]]
disease = "Acute myeloid leukemia"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[26]]
disease = "Acute myeloid leukemia"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[27]]
disease = "Alzheimer's disease"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[28]]
disease = "Type II diabetes mellitus"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[29]]
disease = "Parkinson's disease"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[30]]
disease = "Colorectal cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[31]]
disease = "Chronic myeloid leukemia"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[32]]
disease = "Chronic myeloid leukemia"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[33]]
disease = "Pancreatic cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[34]]
disease = "Colorectal cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

dataset <- resList[[35]]
disease = "Endometrial cancer"
match(disease,rownames(dataset))
rank <- cbind(rank,match(disease,rownames(dataset)))
dataset[disease, 'pG']
pValue <- cbind(pValue,dataset[disease, 'pG'])

boxplot(as.vector(rank),main ="Ranks of target pathways", ylim=c(0, 120))
boxplot(as.vector(pValue),main ="p-values of target pathways", ylim = c(0,1))
hist(rank)
hist(pValue, breaks = 50)
resList$GSE7305
hist(resList$GSE7305)