
path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsPADOG.RData")

load("diseaseList.RData")
dataSets=c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
           "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")

TargetPW = c("Renal cell carcinoma", "Alzheimer's disease", "Thyroid cancer", "Dilated cardiomyopathy",
             "Thyroid cancer", "Colorectal cancer", "Alzheimer's disease", "Alzheimer's disease",
             "Alzheimer's disease", "Prostate cancer", "Prostate cancer", "Colorectal cancer", 
             "Huntington's disease", "Colorectal cancer", "Acute myeloid leukemia", 
             "Renal cell carcinoma", "Pancreatic cancer", "Pancreatic cancer", "Non-small cell lung cancer",
             "Non-small cell lung cancer", "Glioma", "Parkinson's disease", "Parkinson's disease", "Glioma",
             "Acute myeloid leukemia", "Acute myeloid leukemia", "Alzheimer's disease", 
             "Type II diabetes mellitus", "Parkinson's disease", "Colorectal cancer", 
             "Chronic myeloid leukemia", "Chronic myeloid leukemia", "Pancreatic cancer",
             "Colorectal cancer", "Endometrial cancer")


rankWithoutNull <- NULL
pvalueWithoutNull <- NULL
myPADOG <- as.matrix(myPADOGComb)
rownames(myPADOG) <- myPADOG[,'Name']

p <- myPADOG["Renal cell carcinoma",'pPADOG_GSE781']
r <- myPADOG["Renal cell carcinoma",'rPADOG_GSE781']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)


p <- myPADOG["Alzheimer's disease",'pPADOG_GSE1297']
r <- myPADOG["Alzheimer's disease",'rPADOG_GSE1297']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Thyroid cancer",'pPADOG_GSE3467']
r <- myPADOG["Thyroid cancer",'rPADOG_GSE3467']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Dilated cardiomyopathy",'pPADOG_GSE3585']
r <- myPADOG["Dilated cardiomyopathy",'rPADOG_GSE3585']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Thyroid cancer",'pPADOG_GSE3678']
r <- myPADOG["Thyroid cancer",'rPADOG_GSE3678']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Colorectal cancer",'pPADOG_GSE4107']
r <- myPADOG["Colorectal cancer",'rPADOG_GSE4107']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Alzheimer's disease",'pPADOG_GSE5281_EC']
r <- myPADOG["Alzheimer's disease",'rPADOG_GSE5281_EC']
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Alzheimer's disease",'pPADOG_GSE5281_HIP']
r <- myPADOG["Alzheimer's disease",'rPADOG_GSE5281_HIP']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Alzheimer's disease",'pPADOG_GSE5281_VCX']
r <- myPADOG["Alzheimer's disease",'rPADOG_GSE5281_VCX']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Prostate cancer",'pPADOG_GSE6956_AA']
r <- myPADOG["Prostate cancer",'rPADOG_GSE6956_AA']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Prostate cancer",'pPADOG_GSE6956_C']
r <- myPADOG["Prostate cancer",'rPADOG_GSE6956_C']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Colorectal cancer",'pPADOG_GSE8671']
r <- myPADOG["Colorectal cancer",'rPADOG_GSE8671']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Huntington's disease",'pPADOG_GSE8762']
r <- myPADOG["Huntington's disease",'rPADOG_GSE8762']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Colorectal cancer",'pPADOG_GSE9348']
r <- myPADOG["Colorectal cancer",'rPADOG_GSE9348']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Acute myeloid leukemia",'pPADOG_GSE9476']
r <- myPADOG["Acute myeloid leukemia",'rPADOG_GSE9476']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Renal cell carcinoma",'pPADOG_GSE14762']
r <- myPADOG["Renal cell carcinoma",'rPADOG_GSE14762']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Pancreatic cancer",'pPADOG_GSE15471']
r <- myPADOG["Pancreatic cancer",'rPADOG_GSE15471']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Pancreatic cancer",'pPADOG_GSE16515']
r <- myPADOG["Pancreatic cancer",'rPADOG_GSE16515']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Non-small cell lung cancer",'pPADOG_GSE18842']
r <- myPADOG["Non-small cell lung cancer",'rPADOG_GSE18842']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Non-small cell lung cancer",'pPADOG_GSE19188']
r <- myPADOG["Non-small cell lung cancer",'rPADOG_GSE19188']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Glioma",'pPADOG_GSE19728']
r <- myPADOG["Glioma",'rPADOG_GSE19728']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Parkinson's disease",'pPADOG_GSE20153']
r <- myPADOG["Parkinson's disease",'rPADOG_GSE20153']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Parkinson's disease",'pPADOG_GSE20291']
r <- myPADOG["Parkinson's disease",'rPADOG_GSE20291']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Glioma",'pPADOG_GSE21354']
r <- myPADOG["Glioma",'rPADOG_GSE21354']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Acute myeloid leukemia",'pPADOG_GSE14924_CD4']
r <- myPADOG["Acute myeloid leukemia",'rPADOG_GSE14924_CD4']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Acute myeloid leukemia",'pPADOG_GSE14924_CD8']
r <- myPADOG["Acute myeloid leukemia",'rPADOG_GSE14924_CD8']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Alzheimer's disease",'pPADOG_GSE16759']
r <- myPADOG["Alzheimer's disease",'rPADOG_GSE16759']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Type II diabetes mellitus",'pPADOG_GSE19420']
r <- myPADOG["Type II diabetes mellitus",'rPADOG_GSE19420']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Parkinson's disease",'pPADOG_GSE20164']
r <- myPADOG["Parkinson's disease",'rPADOG_GSE20164']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Colorectal cancer",'pPADOG_GSE23878']
r <- myPADOG["Colorectal cancer",'rPADOG_GSE23878']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Chronic myeloid leukemia",'pPADOG_GSE24739_G0']
r <- myPADOG["Chronic myeloid leukemia",'rPADOG_GSE24739_G0']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Chronic myeloid leukemia",'pPADOG_GSE24739_G1']
r <- myPADOG["Chronic myeloid leukemia",'rPADOG_GSE24739_G1']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Pancreatic cancer",'pPADOG_GSE32676']
r <- myPADOG["Pancreatic cancer",'rPADOG_GSE32676']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Colorectal cancer",'pPADOG_GSE4183']
r <- myPADOG["Colorectal cancer",'rPADOG_GSE4183']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- myPADOG["Endometrial cancer",'pPADOG_GSE7305']
r <- myPADOG["Endometrial cancer",'rPADOG_GSE7305']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)


load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsPADOG_corrected.RData")


rankWithNull <- NULL
pvalueWithNull <- NULL
myPADOG_corrected <- as.matrix(myPADOGComb_Corrected)
rownames(myPADOG_corrected) <- myPADOG_corrected[,'Name']

p <- myPADOG_corrected["Renal cell carcinoma",'pPADOG_GSE781']
r <- myPADOG_corrected["Renal cell carcinoma",'rPADOG_GSE781']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)


p <- myPADOG_corrected["Alzheimer's disease",'pPADOG_GSE1297']
r <- myPADOG_corrected["Alzheimer's disease",'rPADOG_GSE1297']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Thyroid cancer",'pPADOG_GSE3467']
r <- myPADOG_corrected["Thyroid cancer",'rPADOG_GSE3467']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Dilated cardiomyopathy",'pPADOG_GSE3585']
r <- myPADOG_corrected["Dilated cardiomyopathy",'rPADOG_GSE3585']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Thyroid cancer",'pPADOG_GSE3678']
r <- myPADOG_corrected["Thyroid cancer",'rPADOG_GSE3678']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Colorectal cancer",'pPADOG_GSE4107']
r <- myPADOG_corrected["Colorectal cancer",'rPADOG_GSE4107']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Alzheimer's disease",'pPADOG_GSE5281_EC']
r <- myPADOG_corrected["Alzheimer's disease",'rPADOG_GSE5281_EC']
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Alzheimer's disease",'pPADOG_GSE5281_HIP']
r <- myPADOG_corrected["Alzheimer's disease",'rPADOG_GSE5281_HIP']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Alzheimer's disease",'pPADOG_GSE5281_VCX']
r <- myPADOG_corrected["Alzheimer's disease",'rPADOG_GSE5281_VCX']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Prostate cancer",'pPADOG_GSE6956_AA']
r <- myPADOG_corrected["Prostate cancer",'rPADOG_GSE6956_AA']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Prostate cancer",'pPADOG_GSE6956_C']
r <- myPADOG_corrected["Prostate cancer",'rPADOG_GSE6956_C']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Colorectal cancer",'pPADOG_GSE8671']
r <- myPADOG_corrected["Colorectal cancer",'rPADOG_GSE8671']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Huntington's disease",'pPADOG_GSE8762']
r <- myPADOG_corrected["Huntington's disease",'rPADOG_GSE8762']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Colorectal cancer",'pPADOG_GSE9348']
r <- myPADOG_corrected["Colorectal cancer",'rPADOG_GSE9348']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Acute myeloid leukemia",'pPADOG_GSE9476']
r <- myPADOG_corrected["Acute myeloid leukemia",'rPADOG_GSE9476']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Renal cell carcinoma",'pPADOG_GSE14762']
r <- myPADOG_corrected["Renal cell carcinoma",'rPADOG_GSE14762']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Pancreatic cancer",'pPADOG_GSE15471']
r <- myPADOG_corrected["Pancreatic cancer",'rPADOG_GSE15471']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Pancreatic cancer",'pPADOG_GSE16515']
r <- myPADOG_corrected["Pancreatic cancer",'rPADOG_GSE16515']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Non-small cell lung cancer",'pPADOG_GSE18842']
r <- myPADOG_corrected["Non-small cell lung cancer",'rPADOG_GSE18842']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Non-small cell lung cancer",'pPADOG_GSE19188']
r <- myPADOG_corrected["Non-small cell lung cancer",'rPADOG_GSE19188']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Glioma",'pPADOG_GSE19728']
r <- myPADOG_corrected["Glioma",'rPADOG_GSE19728']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Parkinson's disease",'pPADOG_GSE20153']
r <- myPADOG_corrected["Parkinson's disease",'rPADOG_GSE20153']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Parkinson's disease",'pPADOG_GSE20291']
r <- myPADOG_corrected["Parkinson's disease",'rPADOG_GSE20291']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Glioma",'pPADOG_GSE21354']
r <- myPADOG_corrected["Glioma",'rPADOG_GSE21354']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Acute myeloid leukemia",'pPADOG_GSE14924_CD4']
r <- myPADOG_corrected["Acute myeloid leukemia",'rPADOG_GSE14924_CD4']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Acute myeloid leukemia",'pPADOG_GSE14924_CD8']
r <- myPADOG_corrected["Acute myeloid leukemia",'rPADOG_GSE14924_CD8']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Alzheimer's disease",'pPADOG_GSE16759']
r <- myPADOG_corrected["Alzheimer's disease",'rPADOG_GSE16759']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Type II diabetes mellitus",'pPADOG_GSE19420']
r <- myPADOG_corrected["Type II diabetes mellitus",'rPADOG_GSE19420']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Parkinson's disease",'pPADOG_GSE20164']
r <- myPADOG_corrected["Parkinson's disease",'rPADOG_GSE20164']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Colorectal cancer",'pPADOG_GSE23878']
r <- myPADOG_corrected["Colorectal cancer",'rPADOG_GSE23878']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Chronic myeloid leukemia",'pPADOG_GSE24739_G0']
r <- myPADOG_corrected["Chronic myeloid leukemia",'rPADOG_GSE24739_G0']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Chronic myeloid leukemia",'pPADOG_GSE24739_G1']
r <- myPADOG_corrected["Chronic myeloid leukemia",'rPADOG_GSE24739_G1']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Pancreatic cancer",'pPADOG_GSE32676']
r <- myPADOG_corrected["Pancreatic cancer",'rPADOG_GSE32676']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Colorectal cancer",'pPADOG_GSE4183']
r <- myPADOG_corrected["Colorectal cancer",'rPADOG_GSE4183']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- myPADOG_corrected["Endometrial cancer",'pPADOG_GSE7305']
r <- myPADOG_corrected["Endometrial cancer",'rPADOG_GSE7305']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

rankWithoutNull <- as.numeric(rankWithoutNull)
pvalueWithoutNull <- as.numeric(pvalueWithoutNull)

rankWithNull <- as.numeric(rankWithNull)
pvalueWithNull <- as.numeric(pvalueWithNull)

median(rankWithoutNull)
median(rankWithNull)
median(pvalueWithoutNull)
median(pvalueWithNull)

names(rankWithoutNull) = TargetPW
names(pvalueWithoutNull)  = TargetPW 
names(rankWithNull) = TargetPW
names(pvalueWithNull) = TargetPW

save(TargetPW, dataSets, rankWithoutNull, pvalueWithoutNull, rankWithNull, pvalueWithNull, file = "PADOGTargetPW.RData")

######################   63 datasets ##############################################

path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("PADOGTargetPW.RData")
rank = rankWithoutNull
pvalue = pvalueWithoutNull

files= c("alzheimerPADOG.RData", "AMLPADOG.RData", "InfluenzaPADOG.RData")
diseases = c("Alzheimer's disease", "Acute myeloid leukemia", "Influenza A")
datasetLength = c(11, 8, 9)
diseaseSymbols = c("path:hsa05010", "path:hsa05221", "path:hsa05164")


for (j in 1:length(files)) {
  load(files[j])
  rankWithoutNull <- rep(NA, datasetLength[j])
  pvalueWithoutNull <- rep(NA, datasetLength[j])
  for (i in 1: datasetLength[j]) {
    print(i)
    temp = PADOGComb[,c(1,i+1)]
    temp = temp[order(temp[,2]),]
    if (!diseaseSymbols[j] %in% rownames(temp)) {
      rankWithoutNull[i] = nrow(temp)+1
      pvalueWithoutNull[i] = 1
    } else{
      rankWithoutNull[i] = which(rownames(temp)== diseaseSymbols[j])
      pvalueWithoutNull[i] = temp[diseaseSymbols[j], 2]
    }
  }
  names(rankWithoutNull) = rep(diseases[j], datasetLength[j])
  names(pvalueWithoutNull) = rep(diseases[j], datasetLength[j])
  rank <- c(rank, rankWithoutNull)
  pvalue <- c(pvalue,pvalueWithoutNull)
}

save(rank, pvalue, file = "63Datasets_PADOG.RData")

######################   54 datasets ##############################################

path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("PADOGTargetPW.RData")
rank = rankWithoutNull
pvalue = pvalueWithoutNull

files= c("alzheimerPADOG.RData", "AMLPADOG.RData")
diseases = c("Alzheimer's disease", "Acute myeloid leukemia")
datasetLength = c(11, 8)
diseaseSymbols = c("path:hsa05010", "path:hsa05221")


for (j in 1:length(files)) {
  load(files[j])
  rankWithoutNull <- rep(NA, datasetLength[j])
  pvalueWithoutNull <- rep(NA, datasetLength[j])
  for (i in 1: datasetLength[j]) {
    print(i)
    temp = PADOGComb[,c(1,i+1)]
    temp = temp[order(temp[,2]),]
    if (!diseaseSymbols[j] %in% rownames(temp)) {
      rankWithoutNull[i] = nrow(temp)+1
      pvalueWithoutNull[i] = 1
    } else{
      rankWithoutNull[i] = which(rownames(temp)== diseaseSymbols[j])
      pvalueWithoutNull[i] = temp[diseaseSymbols[j], 2]
    }
  }
  names(rankWithoutNull) = rep(diseases[j], datasetLength[j])
  names(pvalueWithoutNull) = rep(diseases[j], datasetLength[j])
  rank <- c(rank, rankWithoutNull)
  pvalue <- c(pvalue,pvalueWithoutNull)
}

save(rank, pvalue, file = "54Datasets_PADOG.RData")




















library(vioplot)
vioplot(as.vector(rankWithoutNull), col = "grey")
vioplot(as.vector(pvalueWithoutNull), col = "grey")

vioplot(as.vector(rankWithNull), col = "grey")
vioplot(as.vector(pvalueWithNull), col = "grey")


####################### DISEASE-Wise ####################### 

boxplot(rankWithoutNull[Renal_cell_carcinoma], rankWithNull[Renal_cell_carcinoma],
        main="Rank(s) of Renal cell carcinoma (PADOG)", 
        ylab=paste0("Rank (", length(Renal_cell_carcinoma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Alzheimers_disease], rankWithNull[Alzheimers_disease],
        main="Rank(s) of Alzheimer's disease (PADOG)", 
        ylab=paste0("Rank (", length(Alzheimers_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Thyroid_cancer], rankWithNull[Thyroid_cancer],
        main="Rank(s) of Thyroid cancer (PADOG)", 
        ylab=paste0("Rank (", length(Thyroid_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Dilated_cardiomyopathy], rankWithNull[Dilated_cardiomyopathy],
        main="Rank(s) of Dilated cardiomyopathy (PADOG)", 
        ylab=paste0("Rank (", length(Dilated_cardiomyopathy), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Colorectal_cancer], rankWithNull[Colorectal_cancer],
        main="Rank(s) of Colorectal cancer (PADOG)", 
        ylab=paste0("Rank (", length(Colorectal_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Prostate_cancer], rankWithNull[Prostate_cancer],
        main="Rank(s) of Prostate cancer (PADOG)", 
        ylab=paste0("Rank (", length(Prostate_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Huntingtons_disease], rankWithNull[Huntingtons_disease],
        main="Rank(s) of Huntington's disease (PADOG)", 
        ylab=paste0("Rank (", length(Huntingtons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Acute_myeloid_leukemia], rankWithNull[Acute_myeloid_leukemia],
        main="Rank(s) of Acute myeloid leukemia (PADOG)", 
        ylab=paste0("Rank (", length(Acute_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Pancreatic_cancer], rankWithNull[Pancreatic_cancer],
        main="Rank(s) of Pancreatic cancer (PADOG)", 
        ylab=paste0("Rank (", length(Pancreatic_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Non_small_cell_lung_cancer], rankWithNull[Non_small_cell_lung_cancer],
        main="Rank(s) of Non-small cell lung cancer (PADOG)", 
        ylab=paste0("Rank (", length(Non_small_cell_lung_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Glioma], rankWithNull[Glioma],
        main="Rank(s) of Glioma (PADOG)", 
        ylab=paste0("Rank (", length(Glioma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Parkinsons_disease], rankWithNull[Parkinsons_disease],
        main="Rank(s) of Parkinson's disease (PADOG)", 
        ylab=paste0("Rank (", length(Parkinsons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Type_II_diabetes_mellitus], rankWithNull[Type_II_diabetes_mellitus],
        main="Rank(s) of Type II diabetes mellitus (PADOG)", 
        ylab=paste0("Rank (", length(Type_II_diabetes_mellitus), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Chronic_myeloid_leukemia], rankWithNull[Chronic_myeloid_leukemia],
        main="Rank(s) of Chronic myeloid leukemia (PADOG)", 
        ylab=paste0("Rank (", length(Chronic_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Endometrial_cancer], rankWithNull[Endometrial_cancer],
        main="Rank(s) of Endometrial cancer (PADOG)", 
        ylab=paste0("Rank (", length(Endometrial_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))

boxplot(rankWithoutNull[Cancer], rankWithNull[Cancer],
        main="Rank(s) of Cancer (PADOG)", 
        ylab=paste0("Rank (", length(Cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Non_cancer], rankWithNull[Non_cancer],
        main="Rank(s) of Non-cancer (PADOG)", 
        ylab=paste0("Rank (", length(Non_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull, rankWithNull,
        main="Rank(s) overall (PADOG)", 
        ylab=paste0("Rank (", length(dataSets), " datasets)" ), names=c("Without Danube", "With Danube"))


######### PVALUE ##########

boxplot(pvalueWithoutNull[Renal_cell_carcinoma], pvalueWithNull[Renal_cell_carcinoma],
        main="pValue(s) of Renal cell carcinoma (PADOG)", 
        ylab=paste0("pValue (", length(Renal_cell_carcinoma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Alzheimers_disease], pvalueWithNull[Alzheimers_disease],
        main="pValue(s) of Alzheimer's disease (PADOG)", 
        ylab=paste0("pValue (", length(Alzheimers_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Thyroid_cancer], pvalueWithNull[Thyroid_cancer],
        main="pValue(s) of Thyroid cancer (PADOG)", 
        ylab=paste0("pValue (", length(Thyroid_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Dilated_cardiomyopathy], pvalueWithNull[Dilated_cardiomyopathy],
        main="pValue(s) of Dilated cardiomyopathy (PADOG)", 
        ylab=paste0("pValue (", length(Dilated_cardiomyopathy), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Colorectal_cancer], pvalueWithNull[Colorectal_cancer],
        main="pValue(s) of Colorectal cancer (PADOG)", 
        ylab=paste0("pValue (", length(Colorectal_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Prostate_cancer], pvalueWithNull[Prostate_cancer],
        main="pValue(s) of Prostate cancer (PADOG)", 
        ylab=paste0("pValue (", length(Prostate_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Huntingtons_disease], pvalueWithNull[Huntingtons_disease],
        main="pValue(s) of Huntington's disease (PADOG)", 
        ylab=paste0("pValue (", length(Huntingtons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Acute_myeloid_leukemia], pvalueWithNull[Acute_myeloid_leukemia],
        main="pValue(s) of Acute myeloid leukemia (PADOG)", 
        ylab=paste0("pValue (", length(Acute_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Pancreatic_cancer], pvalueWithNull[Pancreatic_cancer],
        main="pValue(s) of Pancreatic cancer (PADOG)", 
        ylab=paste0("pValue (", length(Pancreatic_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Non_small_cell_lung_cancer], pvalueWithNull[Non_small_cell_lung_cancer],
        main="pValue(s) of Non-small cell lung cancer (PADOG)", 
        ylab=paste0("pValue (", length(Non_small_cell_lung_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Glioma], pvalueWithNull[Glioma],
        main="pValue(s) of Glioma (PADOG)", 
        ylab=paste0("pValue (", length(Glioma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Parkinsons_disease], pvalueWithNull[Parkinsons_disease],
        main="pValue(s) of Parkinson's disease (PADOG)", 
        ylab=paste0("pValue (", length(Parkinsons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Type_II_diabetes_mellitus], pvalueWithNull[Type_II_diabetes_mellitus],
        main="pValue(s) of Type II diabetes mellitus (PADOG)", 
        ylab=paste0("pValue (", length(Type_II_diabetes_mellitus), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Chronic_myeloid_leukemia], pvalueWithNull[Chronic_myeloid_leukemia],
        main="pValue(s) of Chronic myeloid leukemia (PADOG)", 
        ylab=paste0("pValue (", length(Chronic_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Endometrial_cancer], pvalueWithNull[Endometrial_cancer],
        main="pValue(s) of Endometrial cancer (PADOG)", 
        ylab=paste0("pValue (", length(Endometrial_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))

boxplot(pvalueWithoutNull[Cancer], pvalueWithNull[Cancer],
        main="pValue(s) of Cancer (PADOG)", 
        ylab=paste0("pValue (", length(Cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Non_cancer], pvalueWithNull[Non_cancer],
        main="pValue(s) of Non-cancer (PADOG)", 
        ylab=paste0("pValue (", length(Non_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull, pvalueWithNull,
        main="pValue(s) overall (PADOG)", 
        ylab=paste0("pValue (", length(dataSets), " datasets)" ), names=c("Without Danube", "With Danube"))


##################################### BIAS ##################################### 

path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsPADOG.RData")

dataSets=c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
           "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")

TargetPW = c("Renal cell carcinoma", "Alzheimer's disease", "Thyroid cancer", "Dilated cardiomyopathy",
             "Thyroid cancer", "Colorectal cancer", "Alzheimer's disease", "Alzheimer's disease",
             "Alzheimer's disease", "Prostate cancer", "Prostate cancer", "Colorectal cancer", 
             "Huntington's disease", "Colorectal cancer", "Acute myeloid leukemia", 
             "Renal cell carcinoma", "Pancreatic cancer", "Pancreatic cancer", "Non-small cell lung cancer",
             "Non-small cell lung cancer", "Glioma", "Parkinson's disease", "Parkinson's disease", "Glioma",
             "Acute myeloid leukemia", "Acute myeloid leukemia", "Alzheimer's disease", 
             "Type II diabetes mellitus", "Parkinson's disease", "Colorectal cancer", 
             "Chronic myeloid leukemia", "Chronic myeloid leukemia", "Pancreatic cancer",
             "Colorectal cancer", "Endometrial cancer")
load("PADOGTargetPW.RData")

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/PADOG_plus2_15vs15_400_final.RData")

bias0 = vector(mode="numeric", length=nrow(resNullAll))
bias1 = vector(mode="numeric", length=nrow(resNullAll))
names(bias0) = resNullAll$Name
names(bias1) = resNullAll$Name

for (i in 1:nrow(resNullAll)) {
  # print(i)
  x = resNullAll[i,-c(1:3)]
  x =  x[colSums(!is.na(x)) > 0]
  x = as.numeric(x)
  # hist(x, breaks = 50)
  densX = density(x)
  # plot(densX)
  fun = data.frame(x = densX$x, y = densX$y)
  
  area0 = fun[fun$x <= 0.25 & fun$x > 0 ,]
  score0 = (area0$y-1)/((area0$x + 0.1)*(area0$x + 0.1))
  area1 = fun[fun$x >= 0.75 & fun$x < 1 ,]
  score1 = (area1$y-1)/((1.1-area1$x)*(1.1-area1$x))
  
  if(max(score0) >= max(score1)) {
    bias0[i] = max(score0)
    bias1[i] = -1
  }
  if(max(score0) < max(score1)) {
    bias0[i] = -1
    bias1[i] = max(score1)
  }
}


# topBias = c(quantile(bias0)[4], quantile(bias1)[4])
top25bias0 = bias0[which(bias0 >= quantile(bias0)[4] & bias0 > 0)]
top25bias1 = bias1[which(bias1 >= quantile(bias1)[4] & bias1 > 0)]

bias0List = bias0[intersect(names(top25bias0),TargetPW)]
bias1List = bias1[intersect(names(top25bias1),TargetPW)]

rankBias0 = rankWithoutNull[names(rankWithoutNull) %in% names(bias0List)]
rankBias1 = rankWithoutNull[names(rankWithoutNull) %in% names(bias1List)]
pValueBias0 = pvalueWithoutNull[names(rankWithoutNull) %in% names(bias0List)]
pValueBias1 = pvalueWithoutNull[names(rankWithoutNull) %in% names(bias1List)]

rankWithNullBias0 = rankWithNull[names(rankWithoutNull) %in% names(bias0List)]
rankWithNullBias1 = rankWithNull[names(rankWithoutNull) %in% names(bias1List)]
pValueWithNullBias0 = pvalueWithNull[names(rankWithoutNull) %in% names(bias0List)]
pValueWithNullBias1 = pvalueWithNull[names(rankWithoutNull) %in% names(bias1List)]

boxplot(rankBias0, rankWithNullBias0,
        main="Ranks of PW biased towards 0", 
        ylab=paste0("Rank (", length(rankBias0), " datasets)" ), names=c("Without Danube", "With Danube"),
        ylim = c(0,150))
boxplot(rankBias1, rankWithNullBias1,
        main="Ranks of PW biased towards 1", 
        ylab=paste0("Rank (", length(rankBias1), " datasets)" ), names=c("Without Danube", "With Danube"),
        ylim = c(0,150))
boxplot(rankWithoutNull, rankWithNull,
        main="Rank(s) overall", 
        ylab=paste0("Rank (", length(dataSets), " datasets)" ), names=c("Without Danube", "With Danube"),
        ylim = c(0,150))

boxplot(pValueBias0, pValueWithNullBias0,
        main="pValues of PW biased towards 0", 
        ylab=paste0("pValue (", length(pValueBias0), " datasets)" ), names=c("Without Danube", "With Danube"),
        ylim = c(0,1))
boxplot(pValueBias1, pValueWithNullBias1,
        main="pValues of PW biased towards 1", 
        ylab=paste0("pValue (", length(pValueBias1), " datasets)" ), names=c("Without Danube", "With Danube"),
        ylim = c(0,1))
boxplot(pvalueWithoutNull, pvalueWithNull,
        main="pValue(s) overall", 
        ylab=paste0("pValue (", length(dataSets), " datasets)" ), names=c("Without Danube", "With Danube"),
        ylim = c(0,1))

save(TargetPW,dataSets, rankBias0, rankWithNullBias0, rankBias1, rankWithNullBias1, rankWithoutNull, rankWithNull,pValueBias0, 
     pValueWithNullBias0,pValueBias1, pValueWithNullBias1,pvalueWithoutNull, pvalueWithNull, file = "PADOG_bias.RData")


################################# Pearson's Bias ################################# 
library(sfsmisc)
path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("pwNames.RData")

expected <- function(x) {
  y = density(x)
  # plot(y)
  expectedvalue = integrate.xy(y$x,y$x*y$y)
  expectedvalue
}

Pearson <- function(x) {
  # a = expected((x - expected(x))^3)/(sd(x)^3)
  a = expected((x-mean(x))^3)/(sd(x)^3)
  a
} 

load("PADOGTargetPW.RData")

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/PADOG_plus2_15vs15_400_final.RData")

bias = vector(mode="numeric", length=nrow(resNullAll))
names(bias) = resNullAll$Name

for (i in 1:nrow(resNullAll)) {
  x = as.numeric(resNullAll[i,-c(1:3)])
  x =  na.omit(x)
  bias[i] = Pearson(x)
}

bias0List = bias[which(bias >= 0.1)]
bias1List = bias[which(bias <= -0.1)]

bias0List = names(bias0List)
bias1List = names(bias1List)
bias0List %in% pwNames$Name
bias1List %in% pwNames$Name
length(bias0List) + length(bias1List)
nrow(resNullAll)

save(bias0List, bias1List, file = "PADOG_bias.RData")

# bias0List = bias0List[intersect(names(bias0List),TargetPW)]
# bias1List = bias1List[intersect(names(bias1List),TargetPW)]
# 
# rankBias0 = rankWithoutNull[names(rankWithoutNull) %in% names(bias0List)]
# rankBias1 = rankWithoutNull[names(rankWithoutNull) %in% names(bias1List)]
# pValueBias0 = pvalueWithoutNull[names(rankWithoutNull) %in% names(bias0List)]
# pValueBias1 = pvalueWithoutNull[names(rankWithoutNull) %in% names(bias1List)]
# 
# rankWithNullBias0 = rankWithNull[names(rankWithoutNull) %in% names(bias0List)]
# rankWithNullBias1 = rankWithNull[names(rankWithoutNull) %in% names(bias1List)]
# pValueWithNullBias0 = pvalueWithNull[names(rankWithoutNull) %in% names(bias0List)]
# pValueWithNullBias1 = pvalueWithNull[names(rankWithoutNull) %in% names(bias1List)]
# 
# boxplot(rankBias0, rankWithNullBias0,
#         main="Ranks of PW biased towards 0", 
#         ylab=paste0("Rank (", length(rankBias0), " datasets)" ), names=c("Without Danube", "With Danube"),
#         ylim = c(0,150))
# boxplot(rankBias1, rankWithNullBias1,
#         main="Ranks of PW biased towards 1", 
#         ylab=paste0("Rank (", length(rankBias1), " datasets)" ), names=c("Without Danube", "With Danube"),
#         ylim = c(0,150))
# boxplot(rankWithoutNull, rankWithNull,
#         main="Rank(s) overall", 
#         ylab=paste0("Rank (", length(dataSets), " datasets)" ), names=c("Without Danube", "With Danube"),
#         ylim = c(0,150))
# 
# boxplot(pValueBias0, pValueWithNullBias0,
#         main="pValues of PW biased towards 0", 
#         ylab=paste0("pValue (", length(pValueBias0), " datasets)" ), names=c("Without Danube", "With Danube"),
#         ylim = c(0,1))
# boxplot(pValueBias1, pValueWithNullBias1,
#         main="pValues of PW biased towards 1", 
#         ylab=paste0("pValue (", length(pValueBias1), " datasets)" ), names=c("Without Danube", "With Danube"),
#         ylim = c(0,1))
# boxplot(pvalueWithoutNull, pvalueWithNull,
#         main="pValue(s) overall", 
#         ylab=paste0("pValue (", length(dataSets), " datasets)" ), names=c("Without Danube", "With Danube"),
#         ylim = c(0,1))
# 
# 
# 
# save(rankBias0, rankWithNullBias0, rankBias1, rankWithNullBias1, rankWithoutNull, rankWithNull,
#      pValueBias0, pValueWithNullBias0,pValueBias1, pValueWithNullBias1, pvalueWithoutNull, pvalueWithNull, 
#      file = "PADOG_bias.RData")
