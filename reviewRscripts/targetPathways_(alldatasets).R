#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ResultSpia.RData")
#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ResultSpiaNormalDEGenes.RData")
#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/correctIA_ROntoTools_400DE.RData")
path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("diseaseList.RData")
load("alldatasetsIA_400_corrected.RData")
#load("resList_ROntoTools_400DE(tryagain).RData")


rank <- NULL
pValue <- NULL
correctedPvalue <- as.matrix(myIAComb_Corrected)
rownames(correctedPvalue) <- correctedPvalue[,'Name']

p <- correctedPvalue["Renal cell carcinoma",'pIA_GSE781']
r <- correctedPvalue["Renal cell carcinoma",'rIA_GSE781']
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)


p <- correctedPvalue["Alzheimer's disease",'pIA_GSE1297']
r <- correctedPvalue["Alzheimer's disease",'rIA_GSE1297']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Thyroid cancer",'pIA_GSE3467']
r <- correctedPvalue["Thyroid cancer",'rIA_GSE3467']
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Dilated cardiomyopathy",'pIA_GSE3585']
r <- correctedPvalue["Dilated cardiomyopathy",'rIA_GSE3585']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Thyroid cancer",'pIA_GSE3678']
r <- correctedPvalue["Thyroid cancer",'rIA_GSE3678']
dataset <- resList$GSE3678
disease = "Thyroid cancer"
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Colorectal cancer",'pIA_GSE4107']
r <- correctedPvalue["Colorectal cancer",'rIA_GSE4107']
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Alzheimer's disease",'pIA_GSE4107']
r <- correctedPvalue["Alzheimer's disease",'rIA_GSE4107']
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Alzheimer's disease",'pIA_GSE5281_HIP']
r <- correctedPvalue["Alzheimer's disease",'rIA_GSE5281_HIP']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Alzheimer's disease",'pIA_GSE5281_VCX']
r <- correctedPvalue["Alzheimer's disease",'rIA_GSE5281_VCX']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Prostate cancer",'pIA_GSE6956_AA']
r <- correctedPvalue["Prostate cancer",'rIA_GSE6956_AA']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Prostate cancer",'pIA_GSE6956_C']
r <- correctedPvalue["Prostate cancer",'rIA_GSE6956_C']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Colorectal cancer",'pIA_GSE8671']
r <- correctedPvalue["Colorectal cancer",'rIA_GSE8671']
p
r
dataset <- resList$GSE8671
disease = "Colorectal cancer"
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Hungtington's disease",'pIA_GSE8762']
r <- correctedPvalue["Hungtington's disease",'rIA_GSE8762']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Colorectal cancer",'pIA_GSE9348']
r <- correctedPvalue["Colorectal cancer",'rIA_GSE9348']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Acute myeloid leukemia",'pIA_GSE9476']
r <- correctedPvalue["Acute myeloid leukemia",'rIA_GSE9476']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Renal cell carcinoma",'pIA_GSE14762']
r <- correctedPvalue["Renal cell carcinoma",'rIA_GSE14762']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Pancreatic cancer",'pIA_GSE15471']
r <- correctedPvalue["Pancreatic cancer",'rIA_GSE15471']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Pancreatic cancer",'pIA_GSE16515']
r <- correctedPvalue["Pancreatic cancer",'rIA_GSE16515']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Non-small cell lung cancer",'pIA_GSE18842']
r <- correctedPvalue["Non-small cell lung cancer",'rIA_GSE18842']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Non-small cell lung cancer",'pIA_GSE19188']
r <- correctedPvalue["Non-small cell lung cancer",'rIA_GSE19188']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Glioma",'pIA_GSE19728']
r <- correctedPvalue["Glioma",'rIA_GSE19728']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Parkinson's disease",'pIA_GSE20153']
r <- correctedPvalue["Parkinson's disease",'rIA_GSE20153']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Parkinson's disease",'pIA_GSE20291']
r <- correctedPvalue["Parkinson's disease",'rIA_GSE20291']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Glioma",'pIA_GSE21354']
r <- correctedPvalue["Glioma",'rIA_GSE21354']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Acute myeloid leukemia",'pIA_GSE14924_CD4']
r <- correctedPvalue["Acute myeloid leukemia",'rIA_GSE14924_CD4']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Acute myeloid leukemia",'pIA_GSE14924_CD8']
r <- correctedPvalue["Acute myeloid leukemia",'rIA_GSE14924_CD8']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Alzheimer's disease",'pIA_GSE16759']
r <- correctedPvalue["Alzheimer's disease",'rIA_GSE16759']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Type II diabetes mellitus",'pIA_GSE19420']
r <- correctedPvalue["Type II diabetes mellitus",'rIA_GSE19420']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Parkinson's disease",'pIA_GSE20164']
r <- correctedPvalue["Parkinson's disease",'rIA_GSE20164']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Colorectal cancer",'pIA_GSE23878']
r <- correctedPvalue["Colorectal cancer",'rIA_GSE23878']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Chronic myeloid leukemia",'pIA_GSE24739_G0']
r <- correctedPvalue["Chronic myeloid leukemia",'rIA_GSE24739_G0']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Chronic myeloid leukemia",'pIA_GSE24739_G1']
r <- correctedPvalue["Chronic myeloid leukemia",'rIA_GSE24739_G1']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Pancreatic cancer",'pIA_GSE32676']
r <- correctedPvalue["Pancreatic cancer",'rIA_GSE32676']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Colorectal cancer",'pIA_GSE4183']
r <- correctedPvalue["Colorectal cancer",'rIA_GSE4183']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

p <- correctedPvalue["Endometrial cancer",'pIA_GSE7305']
r <- correctedPvalue["Endometrial cancer",'rIA_GSE7305']
p
r
rank <- cbind(rank,r)
pValue <- cbind(pValue,p)

boxplot(as.vector(as.numeric(rank)),main ="Ranks of target pathways", ylim=c(0, 120))
boxplot(as.vector(as.numeric(pValue)),main ="p-values of target pathways", ylim = c(0,1))
hist(rank)
hist(pValue, breaks = 50)
correctIA$GSE7305
hist(correctIA$GSE7305)