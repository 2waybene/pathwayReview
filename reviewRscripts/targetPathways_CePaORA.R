library(vioplot)
path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
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

rankWithoutNull <- rep(NA, 35)
pvalueWithoutNull <- rep(NA, 35)


load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsCePaORA_v2.RData")
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/CePaORA_corrected_v2.RData")

for(i in 1:length(Renal_cell_carcinoma)) {
  temp <- myData[[Renal_cell_carcinoma[i]]]
  rankWithoutNull[Renal_cell_carcinoma[i]] <- which(row.names(temp)=="hsa05211")
  pvalueWithoutNull[Renal_cell_carcinoma[i]] <- temp["hsa05211","min_pValue"]
}

for(i in 1:length(Alzheimers_disease)) {
  temp <- myData[[Alzheimers_disease[i]]]
  rankWithoutNull[Alzheimers_disease[i]] <- which(row.names(temp)=="hsa05010")
  pvalueWithoutNull[Alzheimers_disease[i]] <- temp["hsa05010","min_pValue"]
}


for(i in 1:length(Thyroid_cancer)) {
  temp <- myData[[Thyroid_cancer[i]]]
  rankWithoutNull[Thyroid_cancer[i]] <- which(row.names(temp)=="hsa05216")
  pvalueWithoutNull[Thyroid_cancer[i]] <- temp["hsa05216","min_pValue"]
}

for(i in 1:length(Colorectal_cancer)) {
  temp <- myData[[Colorectal_cancer[i]]]
  rankWithoutNull[Colorectal_cancer[i]] <- which(row.names(temp)=="hsa05210")
  pvalueWithoutNull[Colorectal_cancer[i]] <- temp["hsa05210","min_pValue"]
}

for(i in 1:length(Prostate_cancer)) {
  temp <- myData[[Prostate_cancer[i]]]
  rankWithoutNull[Prostate_cancer[i]] <- which(row.names(temp)=="hsa05215")
  pvalueWithoutNull[Prostate_cancer[i]] <- temp["hsa05215","min_pValue"]
}

for(i in 1:length(Acute_myeloid_leukemia)) {
  temp <- myData[[Acute_myeloid_leukemia[i]]]
  rankWithoutNull[Acute_myeloid_leukemia[i]] <- which(row.names(temp)=="hsa05221")
  pvalueWithoutNull[Acute_myeloid_leukemia[i]] <- temp["hsa05221","min_pValue"]
}

for(i in 1:length(Pancreatic_cancer)) {
  temp <- myData[[Pancreatic_cancer[i]]]
  rankWithoutNull[Pancreatic_cancer[i]] <- which(row.names(temp)=="hsa05212")
  pvalueWithoutNull[Pancreatic_cancer[i]] <- temp["hsa05212","min_pValue"]
}

for(i in 1:length(Non_small_cell_lung_cancer)) {
  temp <- myData[[Non_small_cell_lung_cancer[i]]]
  rankWithoutNull[Non_small_cell_lung_cancer[i]] <- which(row.names(temp)=="hsa05223")
  pvalueWithoutNull[Non_small_cell_lung_cancer[i]] <- temp["hsa05223","min_pValue"]
}

for(i in 1:length(Glioma)) {
  temp <- myData[[Glioma[i]]]
  rankWithoutNull[Glioma[i]] <- which(row.names(temp)=="hsa05214")
  pvalueWithoutNull[Glioma[i]] <- temp["hsa05214","min_pValue"]
}

for(i in 1:length(Parkinsons_disease)) {
  temp <- myData[[Parkinsons_disease[i]]]
  rankWithoutNull[Parkinsons_disease[i]] <- which(row.names(temp)=="hsa05012")
  pvalueWithoutNull[Parkinsons_disease[i]] <- temp["hsa05012","min_pValue"]
}

for(i in 1:length(Chronic_myeloid_leukemia)) {
  temp <- myData[[Chronic_myeloid_leukemia[i]]]
  rankWithoutNull[Chronic_myeloid_leukemia[i]] <-  which(row.names(temp)=="hsa05220")
  pvalueWithoutNull[Chronic_myeloid_leukemia[i]] <- temp["hsa05220","min_pValue"]
}

for(i in 1:length(Type_II_diabetes_mellitus)) {
  temp <- myData[[Type_II_diabetes_mellitus[i]]]
  rankWithoutNull[Type_II_diabetes_mellitus[i]] <- which(row.names(temp)=="hsa04930")
  pvalueWithoutNull[Type_II_diabetes_mellitus[i]] <- temp["hsa04930","min_pValue"]
}


for(i in 1:length(Dilated_cardiomyopathy)) {
  temp <- myData[[Dilated_cardiomyopathy[i]]]
  rankWithoutNull[Dilated_cardiomyopathy[i]] <- which(row.names(temp)=="hsa05414")
  pvalueWithoutNull[Dilated_cardiomyopathy[i]] <- temp["hsa05414","min_pValue"]
}


for(i in 1:length(Huntingtons_disease)) {
  temp <- myData[[Huntingtons_disease[i]]]
  rankWithoutNull[Huntingtons_disease[i]] <- which(row.names(temp)=="hsa05016")
  pvalueWithoutNull[Huntingtons_disease[i]] <- temp["hsa05016","min_pValue"]
}

for(i in 1:length(Endometrial_cancer)) {
  temp <- myData[[Endometrial_cancer[i]]]
  rankWithoutNull[Endometrial_cancer[i]] <- which(row.names(temp)=="hsa05213")
  pvalueWithoutNull[Endometrial_cancer[i]] <- temp["hsa05213","min_pValue"]
}

rankWithNull <- rep(NA, 35)
pvalueWithNull <- rep(NA, 35)

for(i in 1:length(Renal_cell_carcinoma)) {
  rankWithNull[Renal_cell_carcinoma[i]] <- myCePaORAComb_Corrected["hsa05211",36+Renal_cell_carcinoma[i]]
  pvalueWithNull[Renal_cell_carcinoma[i]] <- myCePaORAComb_Corrected["hsa05211",1+Renal_cell_carcinoma[i]]
}

for(i in 1:length(Alzheimers_disease)) {
  rankWithNull[Alzheimers_disease[i]] <- myCePaORAComb_Corrected["hsa05010",36+Alzheimers_disease[i]]
  pvalueWithNull[Alzheimers_disease[i]] <- myCePaORAComb_Corrected["hsa05010",1+Alzheimers_disease[i]]
}


for(i in 1:length(Thyroid_cancer)) {
  
  rankWithNull[Thyroid_cancer[i]] <- myCePaORAComb_Corrected["hsa05216",36+ Thyroid_cancer[i]]
  pvalueWithNull[Thyroid_cancer[i]] <- myCePaORAComb_Corrected["hsa05216",1+Thyroid_cancer[i]]
}

for(i in 1:length(Colorectal_cancer)) {
  #myCePaORAComb_Corrected <- myData[[Colorectal_cancer[i]]]
  rankWithNull[Colorectal_cancer[i]] <- myCePaORAComb_Corrected["hsa05210",36+Colorectal_cancer[i]]
  pvalueWithNull[Colorectal_cancer[i]] <- myCePaORAComb_Corrected["hsa05210",1+Colorectal_cancer[i]]
}

for(i in 1:length(Prostate_cancer)) {
  #myCePaORAComb_Corrected <- myData[[Prostate_cancer[i]]]
  rankWithNull[Prostate_cancer[i]] <- myCePaORAComb_Corrected["hsa05215",36+Prostate_cancer[i]]
  pvalueWithNull[Prostate_cancer[i]] <- myCePaORAComb_Corrected["hsa05215",1+Prostate_cancer[i]]
}

for(i in 1:length(Acute_myeloid_leukemia)) {
  #myCePaORAComb_Corrected <- myData[[Acute_myeloid_leukemia[i]]]
  rankWithNull[Acute_myeloid_leukemia[i]] <- myCePaORAComb_Corrected["hsa05221",36+Acute_myeloid_leukemia[i]]
  pvalueWithNull[Acute_myeloid_leukemia[i]] <- myCePaORAComb_Corrected["hsa05221",1+Acute_myeloid_leukemia[i]]
}

for(i in 1:length(Pancreatic_cancer)) {
  
  rankWithNull[Pancreatic_cancer[i]] <- myCePaORAComb_Corrected["hsa05212",36+Pancreatic_cancer[i]]
  pvalueWithNull[Pancreatic_cancer[i]] <- myCePaORAComb_Corrected["hsa05212",1+Pancreatic_cancer[i]]
}

for(i in 1:length(Non_small_cell_lung_cancer)) {
  
  rankWithNull[Non_small_cell_lung_cancer[i]] <- myCePaORAComb_Corrected["hsa05223",36+Non_small_cell_lung_cancer[i]]
  pvalueWithNull[Non_small_cell_lung_cancer[i]] <- myCePaORAComb_Corrected["hsa05223",1+Non_small_cell_lung_cancer[i]]
}

for(i in 1:length(Glioma)) {
  rankWithNull[Glioma[i]] <- myCePaORAComb_Corrected["hsa05214",36+Glioma[i]]
  pvalueWithNull[Glioma[i]] <- myCePaORAComb_Corrected["hsa05214",1+Glioma[i]]
}

for(i in 1:length(Parkinsons_disease)) {
  rankWithNull[Parkinsons_disease[i]] <- myCePaORAComb_Corrected["hsa05012",36+Parkinsons_disease[i]]
  pvalueWithNull[Parkinsons_disease[i]] <- myCePaORAComb_Corrected["hsa05012",1+Parkinsons_disease[i]]
}

for(i in 1:length(Chronic_myeloid_leukemia)) {
  rankWithNull[Chronic_myeloid_leukemia[i]] <- myCePaORAComb_Corrected["hsa05220",36+Chronic_myeloid_leukemia[i]]
  pvalueWithNull[Chronic_myeloid_leukemia[i]] <- myCePaORAComb_Corrected["hsa05220",Chronic_myeloid_leukemia[i]]
}

for(i in 1:length(Type_II_diabetes_mellitus)) {
  rankWithNull[Type_II_diabetes_mellitus[i]] <- myCePaORAComb_Corrected["hsa04930",36+Type_II_diabetes_mellitus[i]]
  pvalueWithNull[Type_II_diabetes_mellitus[i]] <- myCePaORAComb_Corrected["hsa04930",1+Type_II_diabetes_mellitus[i]]
}


for(i in 1:length(Dilated_cardiomyopathy)) {
  rankWithNull[Dilated_cardiomyopathy[i]] <- myCePaORAComb_Corrected["hsa05414",36+Dilated_cardiomyopathy[i]]
  pvalueWithNull[Dilated_cardiomyopathy[i]] <- myCePaORAComb_Corrected["hsa05414",1+Dilated_cardiomyopathy[i]]
}


for(i in 1:length(Huntingtons_disease)) {
  rankWithNull[Huntingtons_disease[i]] <- myCePaORAComb_Corrected["hsa05016",36+Huntingtons_disease[i]]
  pvalueWithNull[Huntingtons_disease[i]] <- myCePaORAComb_Corrected["hsa05016",1+Huntingtons_disease[i]]
}

for(i in 1:length(Endometrial_cancer)) {
  rankWithNull[Endometrial_cancer[i]] <- myCePaORAComb_Corrected["hsa05213",36+Endometrial_cancer[i]]
  pvalueWithNull[Endometrial_cancer[i]] <- myCePaORAComb_Corrected["hsa05213",1+Endometrial_cancer[i]]
}

names(rankWithoutNull) = TargetPW
names(pvalueWithoutNull)  = TargetPW 
names(rankWithNull) = TargetPW
names(pvalueWithNull) = TargetPW

save(TargetPW, dataSets, rankWithoutNull, pvalueWithoutNull, rankWithNull, pvalueWithNull, file = "CePaORATargetPW.RData")

####################################################################

path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("CePaORATargetPW.RData")
rank = rankWithoutNull
pvalue = pvalueWithoutNull

files= c("alzheimerCePaORA.RData", "AMLCePaORA.RData")
diseases = c("Alzheimer's disease", "Acute myeloid leukemia")
diseaseSymbols = c("hsa05010", "hsa05221")
datasetLength = c(11, 8)

for (j in 1:length(files)) {
  load(files[j])
  rankWithoutNull <- rep(NA, datasetLength[j])
  pvalueWithoutNull <- rep(NA, datasetLength[j])
  for (i in 1: datasetLength[j]) {
    print(i)
    temp = myData[[i]]
    if (!diseaseSymbols[j] %in% rownames(temp)) {
      rankWithoutNull[i] = nrow(temp)+1
      pvalueWithoutNull[i] = 1
    } else{
      rankWithoutNull[i] = temp[diseaseSymbols[j], "rank"]
      pvalueWithoutNull[i] = temp[diseaseSymbols[j], "min_pValue"]
    }
  }
  names(rankWithoutNull) = rep(diseases[j], datasetLength[j])
  names(pvalueWithoutNull) = rep(diseases[j], datasetLength[j])
  rank <- c(rank, rankWithoutNull)
  pvalue <- c(pvalue,pvalueWithoutNull)
}

save(rank, pvalue, file = "54Datasets_CePaORA.RData")
























vioplot(as.vector(as.numeric(rankWithoutNull)), col = "grey", ylim = c(0,120))
title("CePaORA without DANUBE", ylab = "Rank of target pathways")
median(as.vector(as.numeric(rankWithoutNull)))

vioplot(as.vector(as.numeric(pvalueWithoutNull)), col = "grey", ylim = c(0,1))
title("CePaORA without DANUBE", ylab = "pValue of target pathways")
median(as.vector(as.numeric(pvalueWithoutNull)))



vioplot(as.vector(as.numeric(rankWithNull)), col = "grey", ylim = c(0,120))
title("CepaORA with DANUBE", ylab = "Rank of target pathways")
median(as.vector(as.numeric(rankWithNull)))
vioplot(as.vector(as.numeric(pvalueWithNull)), col = "grey", ylim = c(0,1))
title("CePa ORA with DANUBE", ylab = "pValue of target pathways")
median(as.vector(as.numeric(pvalueWithNull)))

col <- c("orange", "blue")
as.vector(rankWithoutNull[Renal_cell_carcinoma])
as.vector(rankIAWithNull[Renal_cell_carcinoma])
plot(as.vector(rankWithoutNull[Renal_cell_carcinoma]), col = col[1], type = 'b', 
     ylim=c(0, 130), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Renal_cell_carcinoma]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Renal_cell_carcinoma]), labels=dataSets[Renal_cell_carcinoma], cex.axis=1)
title("Renal cell carcinoma")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Alzheimers_disease])
as.vector(rankIAWithNull[Alzheimers_disease])
plot(as.vector(rankWithoutNull[Alzheimers_disease]), col = col[1], type = 'b', 
     ylim=c(0, 120), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Alzheimers_disease]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Alzheimers_disease]), labels=dataSets[Alzheimers_disease], cex.axis=1)
title("Alzheimer's disease")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Thyroid_cancer])
as.vector(rankIAWithNull[Thyroid_cancer])
plot(as.vector(rankWithoutNull[Thyroid_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 110), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Thyroid_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Thyroid_cancer]), labels=dataSets[Thyroid_cancer], cex.axis=1)
title("Thyroid cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Dilated_cardiomyopathy])
as.vector(rankIAWithNull[Dilated_cardiomyopathy])
plot(as.vector(rankWithoutNull[Dilated_cardiomyopathy]), col = col[1], type = 'b', 
     ylim=c(0, 90), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Dilated_cardiomyopathy]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Dilated_cardiomyopathy]), labels=dataSets[Dilated_cardiomyopathy], cex.axis=1)
title("Dilated cardiomyopathy")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Colorectal_cancer])
as.vector(rankIAWithNull[Colorectal_cancer])
plot(as.vector(rankWithoutNull[Colorectal_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 90), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Colorectal_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Colorectal_cancer]), labels=dataSets[Colorectal_cancer], cex.axis=1)
title("Colorectal cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Prostate_cancer])
as.vector(rankIAWithNull[Prostate_cancer])
plot(as.vector(rankWithoutNull[Prostate_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 80), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Prostate_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Prostate_cancer]), labels=dataSets[Prostate_cancer], cex.axis=1)
title("Prostate cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Huntingtons_disease])
as.vector(rankIAWithNull[Huntingtons_disease])
plot(as.vector(rankWithoutNull[Huntingtons_disease]), col = col[1], type = 'b', 
     ylim=c(0, 20), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Huntingtons_disease]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Huntingtons_disease]), labels=dataSets[Huntingtons_disease], cex.axis=1)
title("Huntington's disease")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Acute_myeloid_leukemia])
as.vector(rankIAWithNull[Acute_myeloid_leukemia])
plot(as.vector(rankWithoutNull[Acute_myeloid_leukemia]), col = col[1], type = 'b', 
     ylim=c(0, 50), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Acute_myeloid_leukemia]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Acute_myeloid_leukemia]), labels=dataSets[Acute_myeloid_leukemia], cex.axis=1)
title("Acute myeloid leukemia")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Pancreatic_cancer])
as.vector(rankIAWithNull[Pancreatic_cancer])
plot(as.vector(rankWithoutNull[Pancreatic_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 25), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Pancreatic_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Pancreatic_cancer]), labels=dataSets[Pancreatic_cancer], cex.axis=1)
title("Pancreatic cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Non_small_cell_lung_cancer])
as.vector(rankIAWithNull[Non_small_cell_lung_cancer])
plot(as.vector(rankWithoutNull[Non_small_cell_lung_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 50), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Non_small_cell_lung_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Non_small_cell_lung_cancer]), labels=dataSets[Non_small_cell_lung_cancer], cex.axis=1)
title("Non-small cell lung cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Glioma])
as.vector(rankIAWithNull[Glioma])
plot(as.vector(rankWithoutNull[Glioma]), col = col[1], type = 'b', 
     ylim=c(0, 20), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Glioma]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Glioma]), labels=dataSets[Glioma], cex.axis=1)
title("Glioma")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Parkinsons_disease])
as.vector(rankIAWithNull[Parkinsons_disease])
plot(as.vector(rankWithoutNull[Parkinsons_disease]), col = col[1], type = 'b', 
     ylim=c(0, 100), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Parkinsons_disease]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Parkinsons_disease]), labels=dataSets[Parkinsons_disease], cex.axis=1)
title("Parkinson's disease")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Type_II_diabetes_mellitus])
as.vector(rankIAWithNull[Type_II_diabetes_mellitus])
plot(as.vector(rankWithoutNull[Type_II_diabetes_mellitus]), col = col[1], type = 'b', 
     ylim=c(0, 25), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Type_II_diabetes_mellitus]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Type_II_diabetes_mellitus]), labels=dataSets[Type_II_diabetes_mellitus], cex.axis=1)
title("Type II diabetes mellitus")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Chronic_myeloid_leukemia])
as.vector(rankIAWithNull[Chronic_myeloid_leukemia])
plot(as.vector(rankWithoutNull[Chronic_myeloid_leukemia]), col = col[1], type = 'b', 
     ylim=c(0, 95), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Chronic_myeloid_leukemia]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Chronic_myeloid_leukemia]), labels=dataSets[Chronic_myeloid_leukemia], cex.axis=1)
title("Chronic myeloid leukemia")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Endometrial_cancer])
as.vector(rankIAWithNull[Endometrial_cancer])
plot(as.vector(rankWithoutNull[Endometrial_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 110), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankIAWithNull[Endometrial_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Endometrial_cancer]), labels=dataSets[Endometrial_cancer], cex.axis=1)
title("Endometrial cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

dev.off()


################################# BIAS ################################# 
path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
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
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/CePaORA15vs15_400_2000_v2.RData")
load("CePaORATargetPW.RData")

pw = rownames(resNullIA)
pw = gsub("path:","", pw)
resNullCePaORA = resNullCePaORA[intersect(pw,rownames(resNullCePaORA)),]

bias0 = vector(mode="numeric", length=nrow(resNullCePaORA))
bias1 = vector(mode="numeric", length=nrow(resNullCePaORA))
names(bias0) = resNullCePaORA$Name
names(bias1) = resNullCePaORA$Name

for (i in 1:nrow(resNullCePaORA)) {
  print(i)
  x = resNullCePaORA[i,-c(1:3)]
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
     pValueWithNullBias0,pValueBias1, pValueWithNullBias1,pvalueWithoutNull, pvalueWithNull, file = "CePaORA_bias.RData")

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

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ROntoTools15vs15_400_5000_onlyHealthy.RData")
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/CePaORA15vs15_400_2000_v2.RData")
load("CePaORATargetPW.RData")

pw = rownames(resNullIA)
pw = gsub("path:","", pw)
resNullCePaORA = resNullCePaORA[intersect(pw,rownames(resNullCePaORA)),]

bias = vector(mode="numeric", length=nrow(resNullCePaORA))
names(bias) = resNullCePaORA$Name

for (i in 1:nrow(resNullCePaORA)) {
  x = as.numeric(resNullCePaORA[i,-c(1:3)])
  x =  na.omit(x)
  bias[i] = Pearson(x)
}

bias0List = bias[which(bias >= 0.1)]
bias1List = bias[which(bias <= -0.1)]

bias0List = names(bias0List)
bias1List = names(bias1List)
bias0List %in% pwNames$Name
bias1List %in% pwNames$Name

save(bias0List, bias1List, file = "CePaORA_bias.RData")

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
#      file = "CePaORA_bias.RData")


####################### DISEASE-Wise ####################### 

boxplot(rankWithoutNull[Renal_cell_carcinoma], rankWithNull[Renal_cell_carcinoma],
        main="Rank(s) of Renal cell carcinoma (CePaORA)", 
        ylab=paste0("Rank (", length(Renal_cell_carcinoma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Alzheimers_disease], rankWithNull[Alzheimers_disease],
        main="Rank(s) of Alzheimer's disease (CePaORA)", 
        ylab=paste0("Rank (", length(Alzheimers_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Thyroid_cancer], rankWithNull[Thyroid_cancer],
        main="Rank(s) of Thyroid cancer (CePaORA)", 
        ylab=paste0("Rank (", length(Thyroid_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Dilated_cardiomyopathy], rankWithNull[Dilated_cardiomyopathy],
        main="Rank(s) of Dilated cardiomyopathy (CePaORA)", 
        ylab=paste0("Rank (", length(Dilated_cardiomyopathy), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Colorectal_cancer], rankWithNull[Colorectal_cancer],
        main="Rank(s) of Colorectal cancer (CePaORA)", 
        ylab=paste0("Rank (", length(Colorectal_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Prostate_cancer], rankWithNull[Prostate_cancer],
        main="Rank(s) of Prostate cancer (CePaORA)", 
        ylab=paste0("Rank (", length(Prostate_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Huntingtons_disease], rankWithNull[Huntingtons_disease],
        main="Rank(s) of Huntington's disease (CePaORA)", 
        ylab=paste0("Rank (", length(Huntingtons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Acute_myeloid_leukemia], rankWithNull[Acute_myeloid_leukemia],
        main="Rank(s) of Acute myeloid leukemia (CePaORA)", 
        ylab=paste0("Rank (", length(Acute_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Pancreatic_cancer], rankWithNull[Pancreatic_cancer],
        main="Rank(s) of Pancreatic cancer (CePaORA)", 
        ylab=paste0("Rank (", length(Pancreatic_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Non_small_cell_lung_cancer], rankWithNull[Non_small_cell_lung_cancer],
        main="Rank(s) of Non-small cell lung cancer (CePaORA)", 
        ylab=paste0("Rank (", length(Non_small_cell_lung_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Glioma], rankWithNull[Glioma],
        main="Rank(s) of Glioma (CePaORA)", 
        ylab=paste0("Rank (", length(Glioma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Parkinsons_disease], rankWithNull[Parkinsons_disease],
        main="Rank(s) of Parkinson's disease (CePaORA)", 
        ylab=paste0("Rank (", length(Parkinsons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Type_II_diabetes_mellitus], rankWithNull[Type_II_diabetes_mellitus],
        main="Rank(s) of Type II diabetes mellitus (CePaORA)", 
        ylab=paste0("Rank (", length(Type_II_diabetes_mellitus), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Chronic_myeloid_leukemia], rankWithNull[Chronic_myeloid_leukemia],
        main="Rank(s) of Chronic myeloid leukemia (CePaORA)", 
        ylab=paste0("Rank (", length(Chronic_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Endometrial_cancer], rankWithNull[Endometrial_cancer],
        main="Rank(s) of Endometrial cancer (CePaORA)", 
        ylab=paste0("Rank (", length(Endometrial_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))

boxplot(rankWithoutNull[Cancer], rankWithNull[Cancer],
        main="Rank(s) of Cancer (CePaORA)", 
        ylab=paste0("Rank (", length(Cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Non_cancer], rankWithNull[Non_cancer],
        main="Rank(s) of Non-cancer (CePaORA)", 
        ylab=paste0("Rank (", length(Non_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull, rankWithNull,
        main="Rank(s) overall (CePaORA)", 
        ylab=paste0("Rank (", length(dataSets), " datasets)" ), names=c("Without Danube", "With Danube"))


######### PVALUE ##########

boxplot(pvalueWithoutNull[Renal_cell_carcinoma], pvalueWithNull[Renal_cell_carcinoma],
        main="pValue(s) of Renal cell carcinoma (CePaORA)", 
        ylab=paste0("pValue (", length(Renal_cell_carcinoma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Alzheimers_disease], pvalueWithNull[Alzheimers_disease],
        main="pValue(s) of Alzheimer's disease (CePaORA)", 
        ylab=paste0("pValue (", length(Alzheimers_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Thyroid_cancer], pvalueWithNull[Thyroid_cancer],
        main="pValue(s) of Thyroid cancer (CePaORA)", 
        ylab=paste0("pValue (", length(Thyroid_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Dilated_cardiomyopathy], pvalueWithNull[Dilated_cardiomyopathy],
        main="pValue(s) of Dilated cardiomyopathy (CePaORA)", 
        ylab=paste0("pValue (", length(Dilated_cardiomyopathy), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Colorectal_cancer], pvalueWithNull[Colorectal_cancer],
        main="pValue(s) of Colorectal cancer (CePaORA)", 
        ylab=paste0("pValue (", length(Colorectal_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Prostate_cancer], pvalueWithNull[Prostate_cancer],
        main="pValue(s) of Prostate cancer (CePaORA)", 
        ylab=paste0("pValue (", length(Prostate_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Huntingtons_disease], pvalueWithNull[Huntingtons_disease],
        main="pValue(s) of Huntington's disease (CePaORA)", 
        ylab=paste0("pValue (", length(Huntingtons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Acute_myeloid_leukemia], pvalueWithNull[Acute_myeloid_leukemia],
        main="pValue(s) of Acute myeloid leukemia (CePaORA)", 
        ylab=paste0("pValue (", length(Acute_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Pancreatic_cancer], pvalueWithNull[Pancreatic_cancer],
        main="pValue(s) of Pancreatic cancer (CePaORA)", 
        ylab=paste0("pValue (", length(Pancreatic_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Non_small_cell_lung_cancer], pvalueWithNull[Non_small_cell_lung_cancer],
        main="pValue(s) of Non-small cell lung cancer (CePaORA)", 
        ylab=paste0("pValue (", length(Non_small_cell_lung_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Glioma], pvalueWithNull[Glioma],
        main="pValue(s) of Glioma (CePaORA)", 
        ylab=paste0("pValue (", length(Glioma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Parkinsons_disease], pvalueWithNull[Parkinsons_disease],
        main="pValue(s) of Parkinson's disease (CePaORA)", 
        ylab=paste0("pValue (", length(Parkinsons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Type_II_diabetes_mellitus], pvalueWithNull[Type_II_diabetes_mellitus],
        main="pValue(s) of Type II diabetes mellitus (CePaORA)", 
        ylab=paste0("pValue (", length(Type_II_diabetes_mellitus), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Chronic_myeloid_leukemia], pvalueWithNull[Chronic_myeloid_leukemia],
        main="pValue(s) of Chronic myeloid leukemia (CePaORA)", 
        ylab=paste0("pValue (", length(Chronic_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Endometrial_cancer], pvalueWithNull[Endometrial_cancer],
        main="pValue(s) of Endometrial cancer (CePaORA)", 
        ylab=paste0("pValue (", length(Endometrial_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))

boxplot(pvalueWithoutNull[Cancer], pvalueWithNull[Cancer],
        main="pValue(s) of Cancer (CePaORA)", 
        ylab=paste0("pValue (", length(Cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Non_cancer], pvalueWithNull[Non_cancer],
        main="pValue(s) of Non-cancer (CePaORA)", 
        ylab=paste0("pValue (", length(Non_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull, pvalueWithNull,
        main="pValue(s) overall (CePaORA)", 
        ylab=paste0("pValue (", length(dataSets), " datasets)" ), names=c("Without Danube", "With Danube"))

