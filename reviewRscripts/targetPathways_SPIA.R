path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsSPIA.RData")

load("diseaseList.RData")

dataSets=c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
           "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")

names(resList) = dataSets

rankWithoutNull <- rep(NA, 35)
pvalueWithoutNull <- rep(NA, 35)


for(i in 1:length(Renal_cell_carcinoma)) {
  temp <- resList[[Renal_cell_carcinoma[i]]]
  rankWithoutNull[Renal_cell_carcinoma[i]] <- which(row.names(temp) == "Renal cell carcinoma")
  pvalueWithoutNull[Renal_cell_carcinoma[i]] <- temp["Renal cell carcinoma","pG"]
}

for(i in 1:length(Alzheimers_disease)) {
  temp <- resList[[Alzheimers_disease[i]]]
  rankWithoutNull[Alzheimers_disease[i]] <- which(row.names(temp)=="Alzheimer's disease")
  pvalueWithoutNull[Alzheimers_disease[i]] <- temp["Alzheimer's disease","pG"]
}


for(i in 1:length(Thyroid_cancer)) {
  temp <- resList[[Thyroid_cancer[i]]]
  rankWithoutNull[Thyroid_cancer[i]] <- which(row.names(temp)=="Thyroid cancer")
  pvalueWithoutNull[Thyroid_cancer[i]] <- temp["Thyroid cancer","pG"]
}

for(i in 1:length(Colorectal_cancer)) {
  temp <- resList[[Colorectal_cancer[i]]]
  rankWithoutNull[Colorectal_cancer[i]] <- which(row.names(temp)=="Colorectal cancer")
  pvalueWithoutNull[Colorectal_cancer[i]] <- temp["Colorectal cancer","pG"]
}

for(i in 1:length(Prostate_cancer)) {
  temp <- resList[[Prostate_cancer[i]]]
  rankWithoutNull[Prostate_cancer[i]] <- which(row.names(temp)=="Prostate cancer")
  pvalueWithoutNull[Prostate_cancer[i]] <- temp["Prostate cancer","pG"]
}

for(i in 1:length(Acute_myeloid_leukemia)) {
  temp <- resList[[Acute_myeloid_leukemia[i]]]
  rankWithoutNull[Acute_myeloid_leukemia[i]] <- which(row.names(temp)=="Acute myeloid leukemia")
  pvalueWithoutNull[Acute_myeloid_leukemia[i]] <- temp["Acute myeloid leukemia","pG"]
}

for(i in 1:length(Pancreatic_cancer)) {
  temp <- resList[[Pancreatic_cancer[i]]]
  rankWithoutNull[Pancreatic_cancer[i]] <- which(row.names(temp)=="Pancreatic cancer")
  pvalueWithoutNull[Pancreatic_cancer[i]] <- temp["Pancreatic cancer","pG"]
}

for(i in 1:length(Non_small_cell_lung_cancer)) {
  temp <- resList[[Non_small_cell_lung_cancer[i]]]
  rankWithoutNull[Non_small_cell_lung_cancer[i]] <- which(row.names(temp)=="Non-small cell lung cancer")
  pvalueWithoutNull[Non_small_cell_lung_cancer[i]] <- temp["Non-small cell lung cancer","pG"]
}

for(i in 1:length(Glioma)) {
  temp <- resList[[Glioma[i]]]
  rankWithoutNull[Glioma[i]] <- which(row.names(temp)=="Glioma")
  pvalueWithoutNull[Glioma[i]] <- temp["Glioma","pG"]
}

for(i in 1:length(Parkinsons_disease)) {
  temp <- resList[[Parkinsons_disease[i]]]
  rankWithoutNull[Parkinsons_disease[i]] <- which(row.names(temp)=="Parkinson's disease")
  pvalueWithoutNull[Parkinsons_disease[i]] <- temp["Parkinson's disease","pG"]
}

for(i in 1:length(Chronic_myeloid_leukemia)) {
  temp <- resList[[Chronic_myeloid_leukemia[i]]]
  rankWithoutNull[Chronic_myeloid_leukemia[i]] <- which(row.names(temp)=="Chronic myeloid leukemia")
  pvalueWithoutNull[Chronic_myeloid_leukemia[i]] <- temp["Chronic myeloid leukemia","pG"]
}

for(i in 1:length(Type_II_diabetes_mellitus)) {
  temp <- resList[[Type_II_diabetes_mellitus[i]]]
  rankWithoutNull[Type_II_diabetes_mellitus[i]] <- which(row.names(temp)=="Type II diabetes mellitus")
  pvalueWithoutNull[Type_II_diabetes_mellitus[i]] <- temp["Type II diabetes mellitus","pG"]
}


for(i in 1:length(Dilated_cardiomyopathy)) {
  temp <- resList[[Dilated_cardiomyopathy[i]]]
  rankWithoutNull[Dilated_cardiomyopathy[i]] <- which(row.names(temp)=="Dilated cardiomyopathy")
  pvalueWithoutNull[Dilated_cardiomyopathy[i]] <- temp["Dilated cardiomyopathy","pG"]
}


for(i in 1:length(Huntingtons_disease)) {
  temp <- resList[[Huntingtons_disease[i]]]
  rankWithoutNull[Huntingtons_disease[i]] <- which(row.names(temp)=="Huntington's disease")
  pvalueWithoutNull[Huntingtons_disease[i]] <- temp["Huntington's disease","pG"]
}

for(i in 1:length(Endometrial_cancer)) {
  temp <- resList[[Endometrial_cancer[i]]]
  rankWithoutNull[Endometrial_cancer[i]] <- which(row.names(temp)=="Endometrial cancer")
  pvalueWithoutNull[Endometrial_cancer[i]] <- temp["Endometrial cancer","pG"]
}

save(rankWithoutNull, pvalueWithoutNull, file = "SPIATargetPW.RData")

######################   63 datasets ##############################################
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

path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("SPIATargetPW.RData")
names(rankWithoutNull) = TargetPW
names(pvalueWithoutNull) = TargetPW
rank = rankWithoutNull
pvalue = pvalueWithoutNull

files= c("alzheimerSPIA.RData", "AMLSPIA.RData", "InfluenzaSPIA.RData")
diseases = c("Alzheimer's disease", "Acute myeloid leukemia", "Influenza A")
datasetLength = c(11, 8, 9)

for (j in 1:length(files)) {
  load(files[j])
  rankWithoutNull <- rep(NA, datasetLength[j])
  pvalueWithoutNull <- rep(NA, datasetLength[j])
  for (i in 1: datasetLength[j]) {
    print(i)
    temp = resList[[i]]
    if (!diseases[j] %in% rownames(temp)) {
      rankWithoutNull[i] = nrow(temp)+1
      pvalueWithoutNull[i] = 1
    } else{
      rankWithoutNull[i] = which(rownames(temp) == diseases[j])
      pvalueWithoutNull[i] = temp[diseases[j], "pG"]
    }
  }
  names(rankWithoutNull) = rep(diseases[j], datasetLength[j])
  names(pvalueWithoutNull) = rep(diseases[j], datasetLength[j])
  rank <- c(rank, rankWithoutNull)
  pvalue <- c(pvalue,pvalueWithoutNull)
}

save(rank, pvalue, file = "63Datasets_SPIA.RData")

######################   54 datasets ##############################################
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

path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("SPIATargetPW.RData")
names(rankWithoutNull) = TargetPW
names(pvalueWithoutNull) = TargetPW
rank = rankWithoutNull
pvalue = pvalueWithoutNull

files= c("alzheimerSPIA.RData", "AMLSPIA.RData")
diseases = c("Alzheimer's disease", "Acute myeloid leukemia")
datasetLength = c(11, 8)

for (j in 1:length(files)) {
  load(files[j])
  rankWithoutNull <- rep(NA, datasetLength[j])
  pvalueWithoutNull <- rep(NA, datasetLength[j])
  for (i in 1: datasetLength[j]) {
    print(i)
    temp = resList[[i]]
    if (!diseases[j] %in% rownames(temp)) {
      rankWithoutNull[i] = nrow(temp)+1
      pvalueWithoutNull[i] = 1
    } else{
      rankWithoutNull[i] = which(rownames(temp) == diseases[j])
      pvalueWithoutNull[i] = temp[diseases[j], "pG"]
    }
  }
  names(rankWithoutNull) = rep(diseases[j], datasetLength[j])
  names(pvalueWithoutNull) = rep(diseases[j], datasetLength[j])
  rank <- c(rank, rankWithoutNull)
  pvalue <- c(pvalue,pvalueWithoutNull)
}

save(rank, pvalue, file = "54Datasets_SPIA.RData")







