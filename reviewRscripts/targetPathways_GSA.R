#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ResultSpGSA.RData")
#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ResultSpGSANormalDEGenes.RData")
#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/correctIA_GSA_400DE.RData")
path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("diseaseList.RData")
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsGSA.RData")

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
PvalueGSA <- as.matrix(myGSAComb)
rownames(PvalueGSA) <- PvalueGSA[,'Name']

p <- PvalueGSA["Renal cell carcinoma",'pGSA_GSE781']
r <- PvalueGSA["Renal cell carcinoma",'rGSA_GSE781']
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)


p <- PvalueGSA["Alzheimer's disease",'pGSA_GSE1297']
r <- PvalueGSA["Alzheimer's disease",'rGSA_GSE1297']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Thyroid cancer",'pGSA_GSE3467']
r <- PvalueGSA["Thyroid cancer",'rGSA_GSE3467']
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Dilated cardiomyopathy",'pGSA_GSE3585']
r <- PvalueGSA["Dilated cardiomyopathy",'rGSA_GSE3585']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Thyroid cancer",'pGSA_GSE3678']
r <- PvalueGSA["Thyroid cancer",'rGSA_GSE3678']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Colorectal cancer",'pGSA_GSE4107']
r <- PvalueGSA["Colorectal cancer",'rGSA_GSE4107']
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Alzheimer's disease",'pGSA_GSE5281_EC']
r <- PvalueGSA["Alzheimer's disease",'rGSA_GSE5281_EC']
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Alzheimer's disease",'pGSA_GSE5281_HIP']
r <- PvalueGSA["Alzheimer's disease",'rGSA_GSE5281_HIP']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Alzheimer's disease",'pGSA_GSE5281_VCX']
r <- PvalueGSA["Alzheimer's disease",'rGSA_GSE5281_VCX']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Prostate cancer",'pGSA_GSE6956_AA']
r <- PvalueGSA["Prostate cancer",'rGSA_GSE6956_AA']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Prostate cancer",'pGSA_GSE6956_C']
r <- PvalueGSA["Prostate cancer",'rGSA_GSE6956_C']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Colorectal cancer",'pGSA_GSE8671']
r <- PvalueGSA["Colorectal cancer",'rGSA_GSE8671']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Huntington's disease",'pGSA_GSE8762']
r <- PvalueGSA["Huntington's disease",'rGSA_GSE8762']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Colorectal cancer",'pGSA_GSE9348']
r <- PvalueGSA["Colorectal cancer",'rGSA_GSE9348']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Acute myeloid leukemia",'pGSA_GSE9476']
r <- PvalueGSA["Acute myeloid leukemia",'rGSA_GSE9476']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Renal cell carcinoma",'pGSA_GSE14762']
r <- PvalueGSA["Renal cell carcinoma",'rGSA_GSE14762']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Pancreatic cancer",'pGSA_GSE15471']
r <- PvalueGSA["Pancreatic cancer",'rGSA_GSE15471']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Pancreatic cancer",'pGSA_GSE16515']
r <- PvalueGSA["Pancreatic cancer",'rGSA_GSE16515']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Non-small cell lung cancer",'pGSA_GSE18842']
r <- PvalueGSA["Non-small cell lung cancer",'rGSA_GSE18842']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Non-small cell lung cancer",'pGSA_GSE19188']
r <- PvalueGSA["Non-small cell lung cancer",'rGSA_GSE19188']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Glioma",'pGSA_GSE19728']
r <- PvalueGSA["Glioma",'rGSA_GSE19728']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Parkinson's disease",'pGSA_GSE20153']
r <- PvalueGSA["Parkinson's disease",'rGSA_GSE20153']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Parkinson's disease",'pGSA_GSE20291']
r <- PvalueGSA["Parkinson's disease",'rGSA_GSE20291']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Glioma",'pGSA_GSE21354']
r <- PvalueGSA["Glioma",'rGSA_GSE21354']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Acute myeloid leukemia",'pGSA_GSE14924_CD4']
r <- PvalueGSA["Acute myeloid leukemia",'rGSA_GSE14924_CD4']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Acute myeloid leukemia",'pGSA_GSE14924_CD8']
r <- PvalueGSA["Acute myeloid leukemia",'rGSA_GSE14924_CD8']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Alzheimer's disease",'pGSA_GSE16759']
r <- PvalueGSA["Alzheimer's disease",'rGSA_GSE16759']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Type II diabetes mellitus",'pGSA_GSE19420']
r <- PvalueGSA["Type II diabetes mellitus",'rGSA_GSE19420']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Parkinson's disease",'pGSA_GSE20164']
r <- PvalueGSA["Parkinson's disease",'rGSA_GSE20164']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Colorectal cancer",'pGSA_GSE23878']
r <- PvalueGSA["Colorectal cancer",'rGSA_GSE23878']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Chronic myeloid leukemia",'pGSA_GSE24739_G0']
r <- PvalueGSA["Chronic myeloid leukemia",'rGSA_GSE24739_G0']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Chronic myeloid leukemia",'pGSA_GSE24739_G1']
r <- PvalueGSA["Chronic myeloid leukemia",'rGSA_GSE24739_G1']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Pancreatic cancer",'pGSA_GSE32676']
r <- PvalueGSA["Pancreatic cancer",'rGSA_GSE32676']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Colorectal cancer",'pGSA_GSE4183']
r <- PvalueGSA["Colorectal cancer",'rGSA_GSE4183']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

p <- PvalueGSA["Endometrial cancer",'pGSA_GSE7305']
r <- PvalueGSA["Endometrial cancer",'rGSA_GSE7305']
p
r
rankWithoutNull <- cbind(rankWithoutNull,r)
pvalueWithoutNull <- cbind(pvalueWithoutNull,p)

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetGSA_corrected.RData")

rankWithNull <- NULL
pvalueWithNull <- NULL
correctedPvalue <- as.matrix(myGSAComb_Corrected)
rownames(correctedPvalue) <- correctedPvalue[,'Name']

p <- correctedPvalue["Renal cell carcinoma",'pGSA_GSE781']
r <- correctedPvalue["Renal cell carcinoma",'rGSA_GSE781']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)


p <- correctedPvalue["Alzheimer's disease",'pGSA_GSE1297']
r <- correctedPvalue["Alzheimer's disease",'rGSA_GSE1297']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Thyroid cancer",'pGSA_GSE3467']
r <- correctedPvalue["Thyroid cancer",'rGSA_GSE3467']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Dilated cardiomyopathy",'pGSA_GSE3585']
r <- correctedPvalue["Dilated cardiomyopathy",'rGSA_GSE3585']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Thyroid cancer",'pGSA_GSE3678']
r <- correctedPvalue["Thyroid cancer",'rGSA_GSE3678']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Colorectal cancer",'pGSA_GSE4107']
r <- correctedPvalue["Colorectal cancer",'rGSA_GSE4107']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Alzheimer's disease",'pGSA_GSE5281_EC']
r <- correctedPvalue["Alzheimer's disease",'rGSA_GSE5281_EC']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Alzheimer's disease",'pGSA_GSE5281_HIP']
r <- correctedPvalue["Alzheimer's disease",'rGSA_GSE5281_HIP']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Alzheimer's disease",'pGSA_GSE5281_VCX']
r <- correctedPvalue["Alzheimer's disease",'rGSA_GSE5281_VCX']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Prostate cancer",'pGSA_GSE6956_AA']
r <- correctedPvalue["Prostate cancer",'rGSA_GSE6956_AA']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Prostate cancer",'pGSA_GSE6956_C']
r <- correctedPvalue["Prostate cancer",'rGSA_GSE6956_C']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Colorectal cancer",'pGSA_GSE8671']
r <- correctedPvalue["Colorectal cancer",'rGSA_GSE8671']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Huntington's disease",'pGSA_GSE8762']
r <- correctedPvalue["Huntington's disease",'rGSA_GSE8762']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Colorectal cancer",'pGSA_GSE9348']
r <- correctedPvalue["Colorectal cancer",'rGSA_GSE9348']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Acute myeloid leukemia",'pGSA_GSE9476']
r <- correctedPvalue["Acute myeloid leukemia",'rGSA_GSE9476']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Renal cell carcinoma",'pGSA_GSE14762']
r <- correctedPvalue["Renal cell carcinoma",'rGSA_GSE14762']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Pancreatic cancer",'pGSA_GSE15471']
r <- correctedPvalue["Pancreatic cancer",'rGSA_GSE15471']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Pancreatic cancer",'pGSA_GSE16515']
r <- correctedPvalue["Pancreatic cancer",'rGSA_GSE16515']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Non-small cell lung cancer",'pGSA_GSE18842']
r <- correctedPvalue["Non-small cell lung cancer",'rGSA_GSE18842']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Non-small cell lung cancer",'pGSA_GSE19188']
r <- correctedPvalue["Non-small cell lung cancer",'rGSA_GSE19188']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Glioma",'pGSA_GSE19728']
r <- correctedPvalue["Glioma",'rGSA_GSE19728']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Parkinson's disease",'pGSA_GSE20153']
r <- correctedPvalue["Parkinson's disease",'rGSA_GSE20153']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Parkinson's disease",'pGSA_GSE20291']
r <- correctedPvalue["Parkinson's disease",'rGSA_GSE20291']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Glioma",'pGSA_GSE21354']
r <- correctedPvalue["Glioma",'rGSA_GSE21354']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Acute myeloid leukemia",'pGSA_GSE14924_CD4']
r <- correctedPvalue["Acute myeloid leukemia",'rGSA_GSE14924_CD4']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Acute myeloid leukemia",'pGSA_GSE14924_CD8']
r <- correctedPvalue["Acute myeloid leukemia",'rGSA_GSE14924_CD8']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Alzheimer's disease",'pGSA_GSE16759']
r <- correctedPvalue["Alzheimer's disease",'rGSA_GSE16759']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Type II diabetes mellitus",'pGSA_GSE19420']
r <- correctedPvalue["Type II diabetes mellitus",'rGSA_GSE19420']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Parkinson's disease",'pGSA_GSE20164']
r <- correctedPvalue["Parkinson's disease",'rGSA_GSE20164']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Colorectal cancer",'pGSA_GSE23878']
r <- correctedPvalue["Colorectal cancer",'rGSA_GSE23878']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Chronic myeloid leukemia",'pGSA_GSE24739_G0']
r <- correctedPvalue["Chronic myeloid leukemia",'rGSA_GSE24739_G0']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Chronic myeloid leukemia",'pGSA_GSE24739_G1']
r <- correctedPvalue["Chronic myeloid leukemia",'rGSA_GSE24739_G1']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Pancreatic cancer",'pGSA_GSE32676']
r <- correctedPvalue["Pancreatic cancer",'rGSA_GSE32676']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Colorectal cancer",'pGSA_GSE4183']
r <- correctedPvalue["Colorectal cancer",'rGSA_GSE4183']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

p <- correctedPvalue["Endometrial cancer",'pGSA_GSE7305']
r <- correctedPvalue["Endometrial cancer",'rGSA_GSE7305']
p
r
rankWithNull <- cbind(rankWithNull,r)
pvalueWithNull <- cbind(pvalueWithNull,p)

rankWithoutNull = as.numeric(rankWithoutNull) 
pvalueWithoutNull = as.numeric(pvalueWithoutNull) 
rankWithNull = as.numeric(rankWithNull) 
pvalueWithNull = as.numeric(pvalueWithNull) 

names(rankWithoutNull) = TargetPW
names(pvalueWithoutNull)  = TargetPW 
names(rankWithNull) = TargetPW
names(pvalueWithNull) = TargetPW

save(TargetPW, dataSets, rankWithoutNull, pvalueWithoutNull, rankWithNull, pvalueWithNull, file = "GSATargetPW.RData")

######################   63 datasets ##############################################


path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("GSATargetPW.RData")
rank = rankWithoutNull
pvalue = pvalueWithoutNull

files= c("alzheimerGSA.RData", "AMLGSA.RData", "InfluenzaGSA.RData")
diseases = c("Alzheimer's disease", "Acute myeloid leukemia", "Influenza A")
datasetLength = c(11, 8, 9)
diseaseSymbols = c("path:hsa05010", "path:hsa05221", "path:hsa05164")


for (j in 1:length(files)) {
  load(files[j])
  rankWithoutNull <- rep(NA, datasetLength[j])
  pvalueWithoutNull <- rep(NA, datasetLength[j])
  for (i in 1: datasetLength[j]) {
    print(i)
    temp = GSAComb[,c(1,i+1)]
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

save(rank, pvalue, file = "63Datasets_GSA.RData")

######################   54 datasets ##############################################


path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("GSATargetPW.RData")
rank = rankWithoutNull
pvalue = pvalueWithoutNull

files= c("alzheimerGSA.RData", "AMLGSA.RData")
diseases = c("Alzheimer's disease", "Acute myeloid leukemia")
datasetLength = c(11, 8)
diseaseSymbols = c("path:hsa05010", "path:hsa05221")


for (j in 1:length(files)) {
  load(files[j])
  rankWithoutNull <- rep(NA, datasetLength[j])
  pvalueWithoutNull <- rep(NA, datasetLength[j])
  for (i in 1: datasetLength[j]) {
    print(i)
    temp = GSAComb[,c(1,i+1)]
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

save(rank, pvalue, file = "54Datasets_GSA.RData")

























####################### DISEASE-Wise ####################### 

boxplot(rankWithoutNull[Renal_cell_carcinoma], rankWithNull[Renal_cell_carcinoma],
        main="Rank(s) of Renal cell carcinoma (GSA)", 
        ylab=paste0("Rank (", length(Renal_cell_carcinoma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Alzheimers_disease], rankWithNull[Alzheimers_disease],
        main="Rank(s) of Alzheimer's disease (GSA)", 
        ylab=paste0("Rank (", length(Alzheimers_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Thyroid_cancer], rankWithNull[Thyroid_cancer],
        main="Rank(s) of Thyroid cancer (GSA)", 
        ylab=paste0("Rank (", length(Thyroid_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Dilated_cardiomyopathy], rankWithNull[Dilated_cardiomyopathy],
        main="Rank(s) of Dilated cardiomyopathy (GSA)", 
        ylab=paste0("Rank (", length(Dilated_cardiomyopathy), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Colorectal_cancer], rankWithNull[Colorectal_cancer],
        main="Rank(s) of Colorectal cancer (GSA)", 
        ylab=paste0("Rank (", length(Colorectal_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Prostate_cancer], rankWithNull[Prostate_cancer],
        main="Rank(s) of Prostate cancer (GSA)", 
        ylab=paste0("Rank (", length(Prostate_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Huntingtons_disease], rankWithNull[Huntingtons_disease],
        main="Rank(s) of Huntington's disease (GSA)", 
        ylab=paste0("Rank (", length(Huntingtons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Acute_myeloid_leukemia], rankWithNull[Acute_myeloid_leukemia],
        main="Rank(s) of Acute myeloid leukemia (GSA)", 
        ylab=paste0("Rank (", length(Acute_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Pancreatic_cancer], rankWithNull[Pancreatic_cancer],
        main="Rank(s) of Pancreatic cancer (GSA)", 
        ylab=paste0("Rank (", length(Pancreatic_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Non_small_cell_lung_cancer], rankWithNull[Non_small_cell_lung_cancer],
        main="Rank(s) of Non-small cell lung cancer (GSA)", 
        ylab=paste0("Rank (", length(Non_small_cell_lung_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Glioma], rankWithNull[Glioma],
        main="Rank(s) of Glioma (GSA)", 
        ylab=paste0("Rank (", length(Glioma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Parkinsons_disease], rankWithNull[Parkinsons_disease],
        main="Rank(s) of Parkinson's disease (GSA)", 
        ylab=paste0("Rank (", length(Parkinsons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Type_II_diabetes_mellitus], rankWithNull[Type_II_diabetes_mellitus],
        main="Rank(s) of Type II diabetes mellitus (GSA)", 
        ylab=paste0("Rank (", length(Type_II_diabetes_mellitus), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Chronic_myeloid_leukemia], rankWithNull[Chronic_myeloid_leukemia],
        main="Rank(s) of Chronic myeloid leukemia (GSA)", 
        ylab=paste0("Rank (", length(Chronic_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Endometrial_cancer], rankWithNull[Endometrial_cancer],
        main="Rank(s) of Endometrial cancer (GSA)", 
        ylab=paste0("Rank (", length(Endometrial_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))

boxplot(rankWithoutNull[Cancer], rankWithNull[Cancer],
        main="Rank(s) of Cancer (GSA)", 
        ylab=paste0("Rank (", length(Cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Non_cancer], rankWithNull[Non_cancer],
        main="Rank(s) of Non-cancer (GSA)", 
        ylab=paste0("Rank (", length(Non_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull, rankWithNull,
        main="Rank(s) overall (GSA)", 
        ylab=paste0("Rank (", length(dataSets), " datasets)" ), names=c("Without Danube", "With Danube"))


######### PVALUE ##########

boxplot(pvalueWithoutNull[Renal_cell_carcinoma], pvalueWithNull[Renal_cell_carcinoma],
        main="pValue(s) of Renal cell carcinoma (GSA)", 
        ylab=paste0("pValue (", length(Renal_cell_carcinoma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Alzheimers_disease], pvalueWithNull[Alzheimers_disease],
        main="pValue(s) of Alzheimer's disease (GSA)", 
        ylab=paste0("pValue (", length(Alzheimers_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Thyroid_cancer], pvalueWithNull[Thyroid_cancer],
        main="pValue(s) of Thyroid cancer (GSA)", 
        ylab=paste0("pValue (", length(Thyroid_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Dilated_cardiomyopathy], pvalueWithNull[Dilated_cardiomyopathy],
        main="pValue(s) of Dilated cardiomyopathy (GSA)", 
        ylab=paste0("pValue (", length(Dilated_cardiomyopathy), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Colorectal_cancer], pvalueWithNull[Colorectal_cancer],
        main="pValue(s) of Colorectal cancer (GSA)", 
        ylab=paste0("pValue (", length(Colorectal_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Prostate_cancer], pvalueWithNull[Prostate_cancer],
        main="pValue(s) of Prostate cancer (GSA)", 
        ylab=paste0("pValue (", length(Prostate_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Huntingtons_disease], pvalueWithNull[Huntingtons_disease],
        main="pValue(s) of Huntington's disease (GSA)", 
        ylab=paste0("pValue (", length(Huntingtons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Acute_myeloid_leukemia], pvalueWithNull[Acute_myeloid_leukemia],
        main="pValue(s) of Acute myeloid leukemia (GSA)", 
        ylab=paste0("pValue (", length(Acute_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Pancreatic_cancer], pvalueWithNull[Pancreatic_cancer],
        main="pValue(s) of Pancreatic cancer (GSA)", 
        ylab=paste0("pValue (", length(Pancreatic_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Non_small_cell_lung_cancer], pvalueWithNull[Non_small_cell_lung_cancer],
        main="pValue(s) of Non-small cell lung cancer (GSA)", 
        ylab=paste0("pValue (", length(Non_small_cell_lung_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Glioma], pvalueWithNull[Glioma],
        main="pValue(s) of Glioma (GSA)", 
        ylab=paste0("pValue (", length(Glioma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Parkinsons_disease], pvalueWithNull[Parkinsons_disease],
        main="pValue(s) of Parkinson's disease (GSA)", 
        ylab=paste0("pValue (", length(Parkinsons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Type_II_diabetes_mellitus], pvalueWithNull[Type_II_diabetes_mellitus],
        main="pValue(s) of Type II diabetes mellitus (GSA)", 
        ylab=paste0("pValue (", length(Type_II_diabetes_mellitus), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Chronic_myeloid_leukemia], pvalueWithNull[Chronic_myeloid_leukemia],
        main="pValue(s) of Chronic myeloid leukemia (GSA)", 
        ylab=paste0("pValue (", length(Chronic_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Endometrial_cancer], pvalueWithNull[Endometrial_cancer],
        main="pValue(s) of Endometrial cancer (GSA)", 
        ylab=paste0("pValue (", length(Endometrial_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))

boxplot(pvalueWithoutNull[Cancer], pvalueWithNull[Cancer],
        main="pValue(s) of Cancer (GSA)", 
        ylab=paste0("pValue (", length(Cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Non_cancer], pvalueWithNull[Non_cancer],
        main="pValue(s) of Non-cancer (GSA)", 
        ylab=paste0("pValue (", length(Non_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull, pvalueWithNull,
        main="pValue(s) overall (GSA)", 
        ylab=paste0("pValue (", length(dataSets), " datasets)" ), names=c("Without Danube", "With Danube"))

##################################### BIAS ##################################### 
path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsGSA.RData")

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
load("GSATargetPW.RData")

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/GSA_plus2_15vs15_400_final.RData")

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
  plot(densX)
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
     pValueWithNullBias0,pValueBias1, pValueWithNullBias1,pvalueWithoutNull, pvalueWithNull, file = "GSA_bias.RData")

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

load("GSATargetPW.RData")

load("GSA_plus2_15vs15_400_final.RData")
# save(dataSets, TargetPW, rankWithNull, rankWithoutNull, pvalueWithNull, pvalueWithoutNull, file = "ROntoToolsTargetPW.RData")

# dim(resNullAll)

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

save(bias0List, bias1List, file = "GSA_bias.RData")

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
# # library(vioplot)
# # vioplot(rankBias0, col = "grey", c(0,150))
# # vioplot(rankWithNullBias0, col = "grey")
# # vioplot(rankBias1, col = "grey")
# # vioplot(rankWithNullBias1, col = "grey")
# 
# save(rankBias0, rankWithNullBias0, rankBias1, rankWithNullBias1, rankWithoutNull, rankWithNull,
#      pValueBias0, pValueWithNullBias0,pValueBias1, pValueWithNullBias1, pvalueWithoutNull, pvalueWithNull, 
#      file = "GSA_bias.RData")
