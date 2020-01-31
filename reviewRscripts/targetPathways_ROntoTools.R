#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ResultSpia.RData")
#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ResultSpiaNormalDEGenes.RData")
library(vioplot)
path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/resList_ROntoTools_400DE(tryagain).RData")
#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/10_13/ROntoTools_corrected.RData")
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
  pvalueWithoutNull[Renal_cell_carcinoma[i]] <- temp["Renal cell carcinoma","pComb"]
}

for(i in 1:length(Alzheimers_disease)) {
  temp <- resList[[Alzheimers_disease[i]]]
  rankWithoutNull[Alzheimers_disease[i]] <- which(row.names(temp)=="Alzheimer's disease")
  pvalueWithoutNull[Alzheimers_disease[i]] <- temp["Alzheimer's disease","pComb"]
}


for(i in 1:length(Thyroid_cancer)) {
  temp <- resList[[Thyroid_cancer[i]]]
  rankWithoutNull[Thyroid_cancer[i]] <- which(row.names(temp)=="Thyroid cancer")
  pvalueWithoutNull[Thyroid_cancer[i]] <- temp["Thyroid cancer","pComb"]
}

for(i in 1:length(Colorectal_cancer)) {
  temp <- resList[[Colorectal_cancer[i]]]
  rankWithoutNull[Colorectal_cancer[i]] <- which(row.names(temp)=="Colorectal cancer")
  pvalueWithoutNull[Colorectal_cancer[i]] <- temp["Colorectal cancer","pComb"]
}

for(i in 1:length(Prostate_cancer)) {
  temp <- resList[[Prostate_cancer[i]]]
  rankWithoutNull[Prostate_cancer[i]] <- which(row.names(temp)=="Prostate cancer")
  pvalueWithoutNull[Prostate_cancer[i]] <- temp["Prostate cancer","pComb"]
}

for(i in 1:length(Acute_myeloid_leukemia)) {
  temp <- resList[[Acute_myeloid_leukemia[i]]]
  rankWithoutNull[Acute_myeloid_leukemia[i]] <- which(row.names(temp)=="Acute myeloid leukemia")
  pvalueWithoutNull[Acute_myeloid_leukemia[i]] <- temp["Acute myeloid leukemia","pComb"]
}

for(i in 1:length(Pancreatic_cancer)) {
  temp <- resList[[Pancreatic_cancer[i]]]
  rankWithoutNull[Pancreatic_cancer[i]] <- which(row.names(temp)=="Pancreatic cancer")
  pvalueWithoutNull[Pancreatic_cancer[i]] <- temp["Pancreatic cancer","pComb"]
}

for(i in 1:length(Non_small_cell_lung_cancer)) {
  temp <- resList[[Non_small_cell_lung_cancer[i]]]
  rankWithoutNull[Non_small_cell_lung_cancer[i]] <- which(row.names(temp)=="Non-small cell lung cancer")
  pvalueWithoutNull[Non_small_cell_lung_cancer[i]] <- temp["Non-small cell lung cancer","pComb"]
}

for(i in 1:length(Glioma)) {
  temp <- resList[[Glioma[i]]]
  rankWithoutNull[Glioma[i]] <- which(row.names(temp)=="Glioma")
  pvalueWithoutNull[Glioma[i]] <- temp["Glioma","pComb"]
}

for(i in 1:length(Parkinsons_disease)) {
  temp <- resList[[Parkinsons_disease[i]]]
  rankWithoutNull[Parkinsons_disease[i]] <- which(row.names(temp)=="Parkinson's disease")
  pvalueWithoutNull[Parkinsons_disease[i]] <- temp["Parkinson's disease","pComb"]
}

for(i in 1:length(Chronic_myeloid_leukemia)) {
  temp <- resList[[Chronic_myeloid_leukemia[i]]]
  rankWithoutNull[Chronic_myeloid_leukemia[i]] <- which(row.names(temp)=="Chronic myeloid leukemia")
  pvalueWithoutNull[Chronic_myeloid_leukemia[i]] <- temp["Chronic myeloid leukemia","pComb"]
}

for(i in 1:length(Type_II_diabetes_mellitus)) {
  temp <- resList[[Type_II_diabetes_mellitus[i]]]
  rankWithoutNull[Type_II_diabetes_mellitus[i]] <- which(row.names(temp)=="Type II diabetes mellitus")
  pvalueWithoutNull[Type_II_diabetes_mellitus[i]] <- temp["Type II diabetes mellitus","pComb"]
}


for(i in 1:length(Dilated_cardiomyopathy)) {
  temp <- resList[[Dilated_cardiomyopathy[i]]]
  rankWithoutNull[Dilated_cardiomyopathy[i]] <- which(row.names(temp)=="Dilated cardiomyopathy")
  pvalueWithoutNull[Dilated_cardiomyopathy[i]] <- temp["Dilated cardiomyopathy","pComb"]
}


for(i in 1:length(Huntingtons_disease)) {
  temp <- resList[[Huntingtons_disease[i]]]
  rankWithoutNull[Huntingtons_disease[i]] <- which(row.names(temp)=="Huntington's disease")
  pvalueWithoutNull[Huntingtons_disease[i]] <- temp["Huntington's disease","pComb"]
}

for(i in 1:length(Endometrial_cancer)) {
  temp <- resList[[Endometrial_cancer[i]]]
  rankWithoutNull[Endometrial_cancer[i]] <- which(row.names(temp)=="Endometrial cancer")
  pvalueWithoutNull[Endometrial_cancer[i]] <- temp["Endometrial cancer","pComb"]
}


vioplot(rankWithoutNull, col = "grey", ylim = c(0,120))
title("ROntoTools without DANUBE", ylab = "Rank of target pathways")
median(rankWithoutNull)
vioplot(pvalueWithoutNull, col = "grey", ylim = c(0,1))
title("ROntoTools without DANUBE", ylab = "pValue of target pathways")
median(pvalueWithoutNull)


################################CorrectedRontoTools################################################################################################

load("diseaseList.RData")
rankWithNull <- rep(NA, 35)
pvalueWithNull <- rep(NA, 35)

for(i in 1:length(Renal_cell_carcinoma)) {
  rankWithNull[Renal_cell_carcinoma[i]] <- myIAComb_Corrected["Renal cell carcinoma",36+Renal_cell_carcinoma[i]]
  pvalueWithNull[Renal_cell_carcinoma[i]] <- myIAComb_Corrected["Renal cell carcinoma",1+Renal_cell_carcinoma[i]]
}

for(i in 1:length(Alzheimers_disease)) {
  rankWithNull[Alzheimers_disease[i]] <- myIAComb_Corrected["Alzheimer's disease",36+Alzheimers_disease[i]]
  pvalueWithNull[Alzheimers_disease[i]] <- myIAComb_Corrected["Alzheimer's disease",1+Alzheimers_disease[i]]
}


for(i in 1:length(Thyroid_cancer)) {
  
  rankWithNull[Thyroid_cancer[i]] <- myIAComb_Corrected["Thyroid cancer",36+ Thyroid_cancer[i]]
  pvalueWithNull[Thyroid_cancer[i]] <- myIAComb_Corrected["Thyroid cancer",1+Thyroid_cancer[i]]
}

for(i in 1:length(Colorectal_cancer)) {
  #myIAComb_Corrected <- resList[[Colorectal_cancer[i]]]
  rankWithNull[Colorectal_cancer[i]] <- myIAComb_Corrected["Colorectal cancer",36+Colorectal_cancer[i]]
  pvalueWithNull[Colorectal_cancer[i]] <- myIAComb_Corrected["Colorectal cancer",1+Colorectal_cancer[i]]
}

for(i in 1:length(Prostate_cancer)) {
  #myIAComb_Corrected <- resList[[Prostate_cancer[i]]]
  rankWithNull[Prostate_cancer[i]] <- myIAComb_Corrected["Prostate cancer",36+Prostate_cancer[i]]
  pvalueWithNull[Prostate_cancer[i]] <- myIAComb_Corrected["Prostate cancer",1+Prostate_cancer[i]]
}

for(i in 1:length(Acute_myeloid_leukemia)) {
  #myIAComb_Corrected <- resList[[Acute_myeloid_leukemia[i]]]
  rankWithNull[Acute_myeloid_leukemia[i]] <- myIAComb_Corrected["Acute myeloid leukemia",Acute_myeloid_leukemia[i]+36]
  pvalueWithNull[Acute_myeloid_leukemia[i]] <- myIAComb_Corrected["Acute myeloid leukemia",1+Acute_myeloid_leukemia[i]]
}

for(i in 1:length(Pancreatic_cancer)) {
  
  rankWithNull[Pancreatic_cancer[i]] <- myIAComb_Corrected["Pancreatic cancer",36+Pancreatic_cancer[i]]
  pvalueWithNull[Pancreatic_cancer[i]] <- myIAComb_Corrected["Pancreatic cancer",1+Pancreatic_cancer[i]]
}

for(i in 1:length(Non_small_cell_lung_cancer)) {
  
  rankWithNull[Non_small_cell_lung_cancer[i]] <- myIAComb_Corrected["Non-small cell lung cancer",36+Non_small_cell_lung_cancer[i]]
  pvalueWithNull[Non_small_cell_lung_cancer[i]] <- myIAComb_Corrected["Non-small cell lung cancer",1+Non_small_cell_lung_cancer[i]]
}

for(i in 1:length(Glioma)) {
  rankWithNull[Glioma[i]] <- myIAComb_Corrected["Glioma",36+Glioma[i]]
  pvalueWithNull[Glioma[i]] <- myIAComb_Corrected["Glioma",1+Glioma[i]]
}

for(i in 1:length(Parkinsons_disease)) {
  rankWithNull[Parkinsons_disease[i]] <- myIAComb_Corrected["Parkinson's disease",36+Parkinsons_disease[i]]
  pvalueWithNull[Parkinsons_disease[i]] <- myIAComb_Corrected["Parkinson's disease",1+Parkinsons_disease[i]]
}

for(i in 1:length(Chronic_myeloid_leukemia)) {
  rankWithNull[Chronic_myeloid_leukemia[i]] <- myIAComb_Corrected["Chronic myeloid leukemia",36+Chronic_myeloid_leukemia[i]]
  pvalueWithNull[Chronic_myeloid_leukemia[i]] <- myIAComb_Corrected["Chronic myeloid leukemia",1+Chronic_myeloid_leukemia[i]]
}

for(i in 1:length(Type_II_diabetes_mellitus)) {
  rankWithNull[Type_II_diabetes_mellitus[i]] <- myIAComb_Corrected["Type II diabetes mellitus",36+Type_II_diabetes_mellitus[i]]
  pvalueWithNull[Type_II_diabetes_mellitus[i]] <- myIAComb_Corrected["Type II diabetes mellitus",1+Type_II_diabetes_mellitus[i]]
}


for(i in 1:length(Dilated_cardiomyopathy)) {
  rankWithNull[Dilated_cardiomyopathy[i]] <- myIAComb_Corrected["Dilated cardiomyopathy",36+Dilated_cardiomyopathy[i]]
  pvalueWithNull[Dilated_cardiomyopathy[i]] <- myIAComb_Corrected["Dilated cardiomyopathy",1+Dilated_cardiomyopathy[i]]
}


for(i in 1:length(Huntingtons_disease)) {
  rankWithNull[Huntingtons_disease[i]] <- myIAComb_Corrected["Huntington's disease",36+Huntingtons_disease[i]]
  pvalueWithNull[Huntingtons_disease[i]] <- myIAComb_Corrected["Huntington's disease",1+Huntingtons_disease[i]]
}

for(i in 1:length(Endometrial_cancer)) {
  rankWithNull[Endometrial_cancer[i]] <- myIAComb_Corrected["Endometrial cancer",36+Endometrial_cancer[i]]
  pvalueWithNull[Endometrial_cancer[i]] <- myIAComb_Corrected["Endometrial cancer",1+Endometrial_cancer[i]]
}

vioplot(as.vector(as.numeric(rankWithNull)), col = "grey", ylim = c(0,120))
title("R Onto Tools with DANUBE", ylab = "Rank of target pathways")
median(as.vector(as.numeric(rankWithNull)))
vioplot(as.vector(as.numeric(pvalueWithNull)), col = "grey", ylim = c(0,1))
title("R Onto Tools  with DANUBE", ylab = "pValue of target pathways")
median(as.vector(as.numeric(pvalueWithNull)))

names(rankWithoutNull) = TargetPW
names(pvalueWithoutNull)  = TargetPW 
names(rankWithNull) = TargetPW
names(pvalueWithNull) = TargetPW

save(rankWithoutNull, pvalueWithoutNull, rankWithNull, pvalueWithNull, file = "ROntoToolsTargetPW.RData")

################ number of significant pathway ###########

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/resList_ROntoTools_400DE(tryagain).RData")
signPW = vector(mode = "numeric", length = 0)
for (i in 1:length(resList)) {
  signPW = c(signPW, length(which(resList[[i]]$pComb<0.05))/150)
}

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ROntoTools40Datasets.RData")
for (i in 1:length(myData)) {
  signPW = c(signPW, length(which(myData[[i]]$pvalue<0.05))/150)
}

mean(signPW)

######################   63 datasets ##############################################

path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("ROntoToolsTargetPW.RData")
rank = rankWithoutNull
pvalue = pvalueWithoutNull

files= c("alzheimerROntoTools.RData", "AMLROntoTools.RData", "InfluenzaROntoTools.RData")
diseases = c("Alzheimer's disease", "Acute myeloid leukemia", "Influenza A")
datasetLength = c(11, 8, 9)

for (j in 1:length(files)) {
  load(files[j])
  names(rankWithoutNull) = rep(diseases[j], datasetLength[j])
  names(pvalueWithoutNull) = rep(diseases[j], datasetLength[j])
  rank <- c(rank, rankWithoutNull)
  pvalue <- c(pvalue,pvalueWithoutNull)
}

save(rank, pvalue, file = "63Datasets_ROntoTools.RData")

######################   54 datasets ##############################################

path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("ROntoToolsTargetPW.RData")
rank = rankWithoutNull
pvalue = pvalueWithoutNull

files= c("alzheimerROntoTools.RData", "AMLROntoTools.RData")
diseases = c("Alzheimer's disease", "Acute myeloid leukemia")
datasetLength = c(11, 8)

for (j in 1:length(files)) {
  load(files[j])
  names(rankWithoutNull) = rep(diseases[j], datasetLength[j])
  names(pvalueWithoutNull) = rep(diseases[j], datasetLength[j])
  rank <- c(rank, rankWithoutNull)
  pvalue <- c(pvalue,pvalueWithoutNull)
}

save(rank, pvalue, file = "54Datasets_ROntoTools.RData")


















####################### DISEASE-Wise ####################### 

boxplot(rankWithoutNull[Renal_cell_carcinoma], rankWithNull[Renal_cell_carcinoma],
        main="Rank(s) of Renal cell carcinoma (ROntoTools)", 
       ylab=paste0("Rank (", length(Renal_cell_carcinoma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Alzheimers_disease], rankWithNull[Alzheimers_disease],
        main="Rank(s) of Alzheimer's disease (ROntoTools)", 
        ylab=paste0("Rank (", length(Alzheimers_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Thyroid_cancer], rankWithNull[Thyroid_cancer],
        main="Rank(s) of Thyroid cancer (ROntoTools)", 
        ylab=paste0("Rank (", length(Thyroid_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Dilated_cardiomyopathy], rankWithNull[Dilated_cardiomyopathy],
        main="Rank(s) of Dilated cardiomyopathy (ROntoTools)", 
        ylab=paste0("Rank (", length(Dilated_cardiomyopathy), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Colorectal_cancer], rankWithNull[Colorectal_cancer],
        main="Rank(s) of Colorectal cancer (ROntoTools)", 
        ylab=paste0("Rank (", length(Colorectal_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Prostate_cancer], rankWithNull[Prostate_cancer],
        main="Rank(s) of Prostate cancer (ROntoTools)", 
        ylab=paste0("Rank (", length(Prostate_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Huntingtons_disease], rankWithNull[Huntingtons_disease],
        main="Rank(s) of Huntington's disease (ROntoTools)", 
        ylab=paste0("Rank (", length(Huntingtons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Acute_myeloid_leukemia], rankWithNull[Acute_myeloid_leukemia],
        main="Rank(s) of Acute myeloid leukemia (ROntoTools)", 
        ylab=paste0("Rank (", length(Acute_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Pancreatic_cancer], rankWithNull[Pancreatic_cancer],
        main="Rank(s) of Pancreatic cancer (ROntoTools)", 
        ylab=paste0("Rank (", length(Pancreatic_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Non_small_cell_lung_cancer], rankWithNull[Non_small_cell_lung_cancer],
        main="Rank(s) of Non-small cell lung cancer (ROntoTools)", 
        ylab=paste0("Rank (", length(Non_small_cell_lung_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Glioma], rankWithNull[Glioma],
        main="Rank(s) of Glioma (ROntoTools)", 
        ylab=paste0("Rank (", length(Glioma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Parkinsons_disease], rankWithNull[Parkinsons_disease],
        main="Rank(s) of Parkinson's disease (ROntoTools)", 
        ylab=paste0("Rank (", length(Parkinsons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Type_II_diabetes_mellitus], rankWithNull[Type_II_diabetes_mellitus],
        main="Rank(s) of Type II diabetes mellitus (ROntoTools)", 
        ylab=paste0("Rank (", length(Type_II_diabetes_mellitus), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Chronic_myeloid_leukemia], rankWithNull[Chronic_myeloid_leukemia],
        main="Rank(s) of Chronic myeloid leukemia (ROntoTools)", 
        ylab=paste0("Rank (", length(Chronic_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Endometrial_cancer], rankWithNull[Endometrial_cancer],
        main="Rank(s) of Endometrial cancer (ROntoTools)", 
        ylab=paste0("Rank (", length(Endometrial_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))

boxplot(rankWithoutNull[Cancer], rankWithNull[Cancer],
        main="Rank(s) of Cancer (ROntoTools)", 
        ylab=paste0("Rank (", length(Cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Non_cancer], rankWithNull[Non_cancer],
        main="Rank(s) of Non-cancer (ROntoTools)", 
        ylab=paste0("Rank (", length(Non_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull, rankWithNull,
        main="Rank(s) overall (ROntoTools)", 
        ylab=paste0("Rank (", length(dataSets), " datasets)" ), names=c("Without Danube", "With Danube"))


######### PVALUE ##########

boxplot(pvalueWithoutNull[Renal_cell_carcinoma], pvalueWithNull[Renal_cell_carcinoma],
        main="pValue(s) of Renal cell carcinoma (ROntoTools)", 
        ylab=paste0("pValue (", length(Renal_cell_carcinoma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Alzheimers_disease], pvalueWithNull[Alzheimers_disease],
        main="pValue(s) of Alzheimer's disease (ROntoTools)", 
        ylab=paste0("pValue (", length(Alzheimers_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Thyroid_cancer], pvalueWithNull[Thyroid_cancer],
        main="pValue(s) of Thyroid cancer (ROntoTools)", 
        ylab=paste0("pValue (", length(Thyroid_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Dilated_cardiomyopathy], pvalueWithNull[Dilated_cardiomyopathy],
        main="pValue(s) of Dilated cardiomyopathy (ROntoTools)", 
        ylab=paste0("pValue (", length(Dilated_cardiomyopathy), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Colorectal_cancer], pvalueWithNull[Colorectal_cancer],
        main="pValue(s) of Colorectal cancer (ROntoTools)", 
        ylab=paste0("pValue (", length(Colorectal_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Prostate_cancer], pvalueWithNull[Prostate_cancer],
        main="pValue(s) of Prostate cancer (ROntoTools)", 
        ylab=paste0("pValue (", length(Prostate_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Huntingtons_disease], pvalueWithNull[Huntingtons_disease],
        main="pValue(s) of Huntington's disease (ROntoTools)", 
        ylab=paste0("pValue (", length(Huntingtons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Acute_myeloid_leukemia], pvalueWithNull[Acute_myeloid_leukemia],
        main="pValue(s) of Acute myeloid leukemia (ROntoTools)", 
        ylab=paste0("pValue (", length(Acute_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Pancreatic_cancer], pvalueWithNull[Pancreatic_cancer],
        main="pValue(s) of Pancreatic cancer (ROntoTools)", 
        ylab=paste0("pValue (", length(Pancreatic_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Non_small_cell_lung_cancer], pvalueWithNull[Non_small_cell_lung_cancer],
        main="pValue(s) of Non-small cell lung cancer (ROntoTools)", 
        ylab=paste0("pValue (", length(Non_small_cell_lung_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Glioma], pvalueWithNull[Glioma],
        main="pValue(s) of Glioma (ROntoTools)", 
        ylab=paste0("pValue (", length(Glioma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Parkinsons_disease], pvalueWithNull[Parkinsons_disease],
        main="pValue(s) of Parkinson's disease (ROntoTools)", 
        ylab=paste0("pValue (", length(Parkinsons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Type_II_diabetes_mellitus], pvalueWithNull[Type_II_diabetes_mellitus],
        main="pValue(s) of Type II diabetes mellitus (ROntoTools)", 
        ylab=paste0("pValue (", length(Type_II_diabetes_mellitus), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Chronic_myeloid_leukemia], pvalueWithNull[Chronic_myeloid_leukemia],
        main="pValue(s) of Chronic myeloid leukemia (ROntoTools)", 
        ylab=paste0("pValue (", length(Chronic_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Endometrial_cancer], pvalueWithNull[Endometrial_cancer],
        main="pValue(s) of Endometrial cancer (ROntoTools)", 
        ylab=paste0("pValue (", length(Endometrial_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))

boxplot(pvalueWithoutNull[Cancer], pvalueWithNull[Cancer],
        main="pValue(s) of Cancer (ROntoTools)", 
        ylab=paste0("pValue (", length(Cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Non_cancer], pvalueWithNull[Non_cancer],
        main="pValue(s) of Non-cancer (ROntoTools)", 
        ylab=paste0("pValue (", length(Non_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull, pvalueWithNull,
        main="pValue(s) overall (ROntoTools)", 
        ylab=paste0("pValue (", length(dataSets), " datasets)" ), names=c("Without Danube", "With Danube"))


################################# myBIAS ################################# 

path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)


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

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ROntoTools15vs15_400_5000_onlyHealthy.RData")
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ROntoToolsTargetPW.RData")
dim(resNullIA)

bias0 = vector(mode="numeric", length=nrow(resNullIA))
bias1 = vector(mode="numeric", length=nrow(resNullIA))
names(bias0) = resNullIA$Name
names(bias1) = resNullIA$Name

for (i in 1:nrow(resNullIA)) {
  if (i == 41 || i == 42) {} else {
    # print(i)
    x = resNullIA[i,-c(1:3)]
    x =  x[colSums(!is.na(x)) > 0]
    x = as.numeric(x)
    densX = density(x)
    # plot(densX)
    fun = data.frame(x = densX$x, y = densX$y)
  
    # d0 = max(fun[fun$x <= 0.25,2])
    # m0 = fun[which(fun$x <= 0.75 & fun$y == d0),1]
    # d0 = d0-1
    # 
    # d1 = max(fun[fun$x >= 0.75,2])
    # m1 = 1-fun[which(fun$x >= 0.75 & fun$y == d1),1]
    # d1 = d1-1
    # 
    # bias0[i] = d0/(m0*m0)
    # bias1[i] = d1/(m1*m1)
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
}

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

save(rankBias0, rankWithNullBias0, rankBias1, rankWithNullBias1, rankWithoutNull, rankWithNull,pValueBias0, 
     pValueWithNullBias0,pValueBias1, pValueWithNullBias1,pvalueWithoutNull, pvalueWithNull, file = "ROntoTools_bias.RData")

# vioplot(as.vector(rankWithNullBias0), col = "grey", ylim=c(0, 120))
# title("Ranks of target pathways (biased towards 0)")
# median(rankWithNullBias0)
# vioplot(as.vector(pValueWithNullBias0), col = "grey",ylim=c(0, 1))
# title("pValue of target pathways (biased towards 0)")
# median(pValueWithNullBias0)
# 
# 
# vioplot(as.vector(rankWithNullBias1), col = "grey",ylim=c(0, 120))
# title("Ranks of target pathways (biased towards 1)")
# median(rankWithNullBias1)
# vioplot(as.vector(pValueWithNullBias1), col = "grey", ylim=c(0, 1))
# title("pValue of target pathways (biased towards 1)")
# median(pValueWithNullBias1)


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
  a = expected((x-mean(x))^3)/(sd(x)^3)
  a
} 

load("ROntoTools15vs15_400_5000_onlyHealthy.RData")
load("ROntoToolsTargetPW.RData")

bias = vector(mode="numeric", length=nrow(resNullIA))
names(bias) = resNullIA$Name

for (i in 1:nrow(resNullIA)) {
  if (i == 41 || i == 42) {} else {
    x = as.numeric(resNullIA[i,-c(1:3)])
    x =  na.omit(x)
    bias[i] = Pearson(x)
  }
}


bias0List = bias[which(bias >= 0.1)]
bias1List = bias[which(bias <= -0.1)]

bias0List = names(bias0List)
bias1List = names(bias1List)

save(bias0List, bias1List, file = "ROntoTools_bias.RData")


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

# library(vioplot)
# vioplot(rankBias0, col = "grey", c(0,150))
# vioplot(rankWithNullBias0, col = "grey")
# vioplot(rankBias1, col = "grey")
# vioplot(rankWithNullBias1, col = "grey")

# save(rankBias0, rankWithNullBias0, rankBias1, rankWithNullBias1, rankWithoutNull, rankWithNull,
#      pValueBias0, pValueWithNullBias0,pValueBias1, pValueWithNullBias1, pvalueWithoutNull, pvalueWithNull, 
#      file = "ROntoTools_bias.RData")






