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

rankPathNetWithoutNull <- rep(NA, 35)
pValuePathNetWithoutNull <- rep(NA, 35)


load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsPathNet (3rd run).RData")

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
  rankPathNetWithoutNull[Chronic_myeloid_leukemia[i]] <- rownames(myData[[Chronic_myeloid_leukemia[i]]][myData[[Chronic_myeloid_leukemia[i]]][,1]=="Chronic myeloid leukemia",])
  pValuePathNetWithoutNull[Chronic_myeloid_leukemia[i]] <- myData[[Chronic_myeloid_leukemia[i]]][myData[[Chronic_myeloid_leukemia[i]]][,1]=="Chronic myeloid leukemia","p_PathNet"]
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

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/PathNet_corrected (3rd run).RData")

rankPathNetWithNull <- rep(NA, 35)
pValuePathNetWithNull <- rep(NA, 35)

for(i in 1:length(Renal_cell_carcinoma)) {
  rankPathNetWithNull[Renal_cell_carcinoma[i]] <- myPathNetComb_Corrected["Renal cell carcinoma",36+Renal_cell_carcinoma[i]]
  pValuePathNetWithNull[Renal_cell_carcinoma[i]] <- myPathNetComb_Corrected["Renal cell carcinoma",1+Renal_cell_carcinoma[i]]
}

for(i in 1:length(Alzheimers_disease)) {
  rankPathNetWithNull[Alzheimers_disease[i]] <- myPathNetComb_Corrected["Alzheimers disease",36+Alzheimers_disease[i]]
  pValuePathNetWithNull[Alzheimers_disease[i]] <- myPathNetComb_Corrected["Alzheimers disease",1+Alzheimers_disease[i]]
}


for(i in 1:length(Thyroid_cancer)) {
  
  rankPathNetWithNull[Thyroid_cancer[i]] <- myPathNetComb_Corrected["Thyroid cancer",36+ Thyroid_cancer[i]]
  pValuePathNetWithNull[Thyroid_cancer[i]] <- myPathNetComb_Corrected["Thyroid cancer",1+Thyroid_cancer[i]]
}

for(i in 1:length(Colorectal_cancer)) {
  #myPathNetComb_Corrected <- myData[[Colorectal_cancer[i]]]
  rankPathNetWithNull[Colorectal_cancer[i]] <- myPathNetComb_Corrected["Colorectal cancer",36+Colorectal_cancer[i]]
  pValuePathNetWithNull[Colorectal_cancer[i]] <- myPathNetComb_Corrected["Colorectal cancer",1+Colorectal_cancer[i]]
}

for(i in 1:length(Prostate_cancer)) {
  #myPathNetComb_Corrected <- myData[[Prostate_cancer[i]]]
  rankPathNetWithNull[Prostate_cancer[i]] <- myPathNetComb_Corrected["Prostate cancer",36+Prostate_cancer[i]]
  pValuePathNetWithNull[Prostate_cancer[i]] <- myPathNetComb_Corrected["Prostate cancer",1+Prostate_cancer[i]]
}

for(i in 1:length(Acute_myeloid_leukemia)) {
  #myPathNetComb_Corrected <- myData[[Acute_myeloid_leukemia[i]]]
  rankPathNetWithNull[Acute_myeloid_leukemia[i]] <- myPathNetComb_Corrected["Acute myeloid leukemia",Acute_myeloid_leukemia[i]+36]
  pValuePathNetWithNull[Acute_myeloid_leukemia[i]] <- myPathNetComb_Corrected["Acute myeloid leukemia",1+Acute_myeloid_leukemia[i]]
}

for(i in 1:length(Pancreatic_cancer)) {
  
  rankPathNetWithNull[Pancreatic_cancer[i]] <- myPathNetComb_Corrected["Pancreatic cancer",36+Pancreatic_cancer[i]]
  pValuePathNetWithNull[Pancreatic_cancer[i]] <- myPathNetComb_Corrected["Pancreatic cancer",1+Pancreatic_cancer[i]]
}

for(i in 1:length(Non_small_cell_lung_cancer)) {
  
  rankPathNetWithNull[Non_small_cell_lung_cancer[i]] <- myPathNetComb_Corrected["Non-small cell lung cancer",36+Non_small_cell_lung_cancer[i]]
  pValuePathNetWithNull[Non_small_cell_lung_cancer[i]] <- myPathNetComb_Corrected["Non-small cell lung cancer",1+Non_small_cell_lung_cancer[i]]
}

for(i in 1:length(Glioma)) {
  rankPathNetWithNull[Glioma[i]] <- myPathNetComb_Corrected["Glioma",36+Glioma[i]]
  pValuePathNetWithNull[Glioma[i]] <- myPathNetComb_Corrected["Glioma",1+Glioma[i]]
}

for(i in 1:length(Parkinsons_disease)) {
  rankPathNetWithNull[Parkinsons_disease[i]] <- myPathNetComb_Corrected["Parkinsons disease",36+Parkinsons_disease[i]]
  pValuePathNetWithNull[Parkinsons_disease[i]] <- myPathNetComb_Corrected["Parkinsons disease",1+Parkinsons_disease[i]]
}

for(i in 1:length(Chronic_myeloid_leukemia)) {
  rankPathNetWithNull[Chronic_myeloid_leukemia[i]] <- myPathNetComb_Corrected["Chronic myeloid leukemia",36+Chronic_myeloid_leukemia[i]]
  pValuePathNetWithNull[Chronic_myeloid_leukemia[i]] <- myPathNetComb_Corrected["Chronic myeloid leukemia",1+Chronic_myeloid_leukemia[i]]
}

for(i in 1:length(Type_II_diabetes_mellitus)) {
  rankPathNetWithNull[Type_II_diabetes_mellitus[i]] <- myPathNetComb_Corrected["Type II diabetes mellitus",36+Type_II_diabetes_mellitus[i]]
  pValuePathNetWithNull[Type_II_diabetes_mellitus[i]] <- myPathNetComb_Corrected["Type II diabetes mellitus",1+Type_II_diabetes_mellitus[i]]
}


for(i in 1:length(Dilated_cardiomyopathy)) {
  rankPathNetWithNull[Dilated_cardiomyopathy[i]] <- myPathNetComb_Corrected["Dilated cardiomyopathy",36+Dilated_cardiomyopathy[i]]
  pValuePathNetWithNull[Dilated_cardiomyopathy[i]] <- myPathNetComb_Corrected["Dilated cardiomyopathy",1+Dilated_cardiomyopathy[i]]
}


for(i in 1:length(Huntingtons_disease)) {
  rankPathNetWithNull[Huntingtons_disease[i]] <- myPathNetComb_Corrected["Huntingtons disease",36+Huntingtons_disease[i]]
  pValuePathNetWithNull[Huntingtons_disease[i]] <- myPathNetComb_Corrected["Huntingtons disease",1+Huntingtons_disease[i]]
}

for(i in 1:length(Endometrial_cancer)) {
  rankPathNetWithNull[Endometrial_cancer[i]] <- myPathNetComb_Corrected["Endometrial cancer",36+Endometrial_cancer[i]]
  pValuePathNetWithNull[Endometrial_cancer[i]] <- myPathNetComb_Corrected["Endometrial cancer",1+Endometrial_cancer[i]]
}

rankWithoutNull = as.numeric(rankPathNetWithoutNull)
pvalueWithoutNull = as.numeric(pValuePathNetWithoutNull)
rankWithNull = as.numeric(rankPathNetWithNull)
pvalueWithNull = as.numeric(pValuePathNetWithNull)

names(rankWithoutNull) = TargetPW
names(pvalueWithoutNull)  = TargetPW 
names(rankWithNull) = TargetPW
names(pvalueWithNull) = TargetPW

save(rankWithoutNull, pvalueWithoutNull, rankWithNull, pvalueWithNull, file = "PathNetTargetPW.RData")

###########################    63 datasets   #########################################

path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("PathNetTargetPW.RData")
rank = rankWithoutNull
pvalue = pvalueWithoutNull

files= c("alzheimerPathNet.RData", "AMLPathNet.RData")
diseases = c("Alzheimers disease", "Acute myeloid leukemia")
datasetLength = c(11, 8)

for (j in 1:length(files)) {
  load(files[j])
  rankWithoutNull <- rep(NA, datasetLength[j])
  pvalueWithoutNull <- rep(NA, datasetLength[j])
  for (i in 1: datasetLength[j]) {
    print(i)
    temp = resList[[i]]
    rownames(temp) = temp$Name
    if (!diseases[j] %in% rownames(temp)) {
      rankWithoutNull[i] = nrow(temp)+1
      pvalueWithoutNull[i] = 1
    } else{
      rankWithoutNull[i] = which(rownames(temp) == diseases[j])
      pvalueWithoutNull[i] = temp[diseases[j], "p_PathNet"]
    }
  }
  names(rankWithoutNull) = rep(diseases[j], datasetLength[j])
  names(pvalueWithoutNull) = rep(diseases[j], datasetLength[j])
  rank <- c(rank, rankWithoutNull)
  pvalue <- c(pvalue,pvalueWithoutNull)
}

save(rank, pvalue, file = "63Datasets_PathNet.RData")





























boxplot(as.vector(as.numeric(rankPathNetWithNull)),main ="Ranks of target pathways with Danube", ylim=c(0, 120))
boxplot(as.vector(as.numeric(rankPathNetWithoutNull)),main ="Ranks of target pathways without Danube", ylim=c(0, 120))
boxplot(as.vector(as.numeric(pValuePathNetWithNull)),main ="p-values of target pathways", ylim = c(0,1))
boxplot(as.vector(as.numeric(pValuePathNetWithoutNull)),main ="p-values of target pathways", ylim = c(0,1))


vioplot(as.vector(as.numeric(rankPathNetWithNull)), col = "grey", ylim = c(0,120))
title("PathNet with DANUBE", ylab = "Rank of target pathways")
median(as.numeric(rankPathNetWithNull))
vioplot(as.vector(as.numeric(pValuePathNetWithNull)), col = "grey", ylim = c(0,1))
title("PathNet with DANUBE", ylab = "pValue of target pathways")
median(as.vector(as.numeric(pValuePathNetWithNull)))


vioplot(as.vector(as.numeric(rankPathNetWithoutNull)), col = "grey", ylim = c(0,120))
title("PathNet without DANUBE", ylab = "Rank of target pathways")
median(as.vector(as.numeric(rankPathNetWithoutNull)))
vioplot(as.vector(as.numeric(pValuePathNetWithoutNull)), col = "grey", ylim = c(0,1))
title("PathNet without DANUBE", ylab = "pValie of target pathways")
median(as.vector(as.numeric(pValuePathNetWithoutNull)))

col <- c("orange", "blue")

as.vector(rankPathNetWithoutNull[Renal_cell_carcinoma])
as.vector(rankPathNetWithNull[Renal_cell_carcinoma])
plot(as.vector(rankPathNetWithoutNull[Renal_cell_carcinoma]), col = col[1], type = 'b', 
     ylim=c(0, 130), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Renal_cell_carcinoma]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Renal_cell_carcinoma]), labels=dataSets[Renal_cell_carcinoma], cex.axis=1)
title("Renal cell carcinoma")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankPathNetWithoutNull[Alzheimers_disease])
as.vector(rankPathNetWithNull[Alzheimers_disease])
plot(as.vector(rankPathNetWithoutNull[Alzheimers_disease]), col = col[1], type = 'b', 
     ylim=c(0, 120), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Alzheimers_disease]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Alzheimers_disease]), labels=dataSets[Alzheimers_disease], cex.axis=1)
title("Alzheimer's disease")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankPathNetWithoutNull[Thyroid_cancer])
as.vector(rankPathNetWithNull[Thyroid_cancer])
plot(as.vector(rankPathNetWithoutNull[Thyroid_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 110), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Thyroid_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Thyroid_cancer]), labels=dataSets[Thyroid_cancer], cex.axis=1)
title("Thyroid cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankPathNetWithoutNull[Dilated_cardiomyopathy])
as.vector(rankPathNetWithNull[Dilated_cardiomyopathy])
plot(as.vector(rankPathNetWithoutNull[Dilated_cardiomyopathy]), col = col[1], type = 'b', 
     ylim=c(0, 90), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Dilated_cardiomyopathy]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Dilated_cardiomyopathy]), labels=dataSets[Dilated_cardiomyopathy], cex.axis=1)
title("Dilated cardiomyopathy")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankPathNetWithoutNull[Colorectal_cancer])
as.vector(rankPathNetWithNull[Colorectal_cancer])
plot(as.vector(rankPathNetWithoutNull[Colorectal_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 90), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Colorectal_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Colorectal_cancer]), labels=dataSets[Colorectal_cancer], cex.axis=1)
title("Colorectal cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankPathNetWithoutNull[Prostate_cancer])
as.vector(rankPathNetWithNull[Prostate_cancer])
plot(as.vector(rankPathNetWithoutNull[Prostate_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 80), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Prostate_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Prostate_cancer]), labels=dataSets[Prostate_cancer], cex.axis=1)
title("Prostate cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankPathNetWithoutNull[Huntingtons_disease])
as.vector(rankPathNetWithNull[Huntingtons_disease])
plot(as.vector(rankPathNetWithoutNull[Huntingtons_disease]), col = col[1], type = 'b', 
     ylim=c(0, 20), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Huntingtons_disease]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Huntingtons_disease]), labels=dataSets[Huntingtons_disease], cex.axis=1)
title("Huntington's disease")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankPathNetWithoutNull[Acute_myeloid_leukemia])
as.vector(rankPathNetWithNull[Acute_myeloid_leukemia])
plot(as.vector(rankPathNetWithoutNull[Acute_myeloid_leukemia]), col = col[1], type = 'b', 
     ylim=c(0, 50), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Acute_myeloid_leukemia]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Acute_myeloid_leukemia]), labels=dataSets[Acute_myeloid_leukemia], cex.axis=1)
title("Acute myeloid leukemia")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankPathNetWithoutNull[Pancreatic_cancer])
as.vector(rankPathNetWithNull[Pancreatic_cancer])
plot(as.vector(rankPathNetWithoutNull[Pancreatic_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 25), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Pancreatic_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Pancreatic_cancer]), labels=dataSets[Pancreatic_cancer], cex.axis=1)
title("Pancreatic cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankPathNetWithoutNull[Non_small_cell_lung_cancer])
as.vector(rankPathNetWithNull[Non_small_cell_lung_cancer])
plot(as.vector(rankPathNetWithoutNull[Non_small_cell_lung_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 50), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Non_small_cell_lung_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Non_small_cell_lung_cancer]), labels=dataSets[Non_small_cell_lung_cancer], cex.axis=1)
title("Non-small cell lung cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankPathNetWithoutNull[Glioma])
as.vector(rankPathNetWithNull[Glioma])
plot(as.vector(rankPathNetWithoutNull[Glioma]), col = col[1], type = 'b', 
     ylim=c(0, 20), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Glioma]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Glioma]), labels=dataSets[Glioma], cex.axis=1)
title("Glioma")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankPathNetWithoutNull[Parkinsons_disease])
as.vector(rankPathNetWithNull[Parkinsons_disease])
plot(as.vector(rankPathNetWithoutNull[Parkinsons_disease]), col = col[1], type = 'b', 
     ylim=c(0, 100), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Parkinsons_disease]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Parkinsons_disease]), labels=dataSets[Parkinsons_disease], cex.axis=1)
title("Parkinson's disease")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankPathNetWithoutNull[Type_II_diabetes_mellitus])
as.vector(rankPathNetWithNull[Type_II_diabetes_mellitus])
plot(as.vector(rankPathNetWithoutNull[Type_II_diabetes_mellitus]), col = col[1], type = 'b', 
     ylim=c(0, 25), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Type_II_diabetes_mellitus]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Type_II_diabetes_mellitus]), labels=dataSets[Type_II_diabetes_mellitus], cex.axis=1)
title("Type II diabetes mellitus")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankPathNetWithoutNull[Chronic_myeloid_leukemia])
as.vector(rankPathNetWithNull[Chronic_myeloid_leukemia])
plot(as.vector(rankPathNetWithoutNull[Chronic_myeloid_leukemia]), col = col[1], type = 'b', 
     ylim=c(0, 95), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Chronic_myeloid_leukemia]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Chronic_myeloid_leukemia]), labels=dataSets[Chronic_myeloid_leukemia], cex.axis=1)
title("Chronic myeloid leukemia")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankPathNetWithoutNull[Endometrial_cancer])
as.vector(rankPathNetWithNull[Endometrial_cancer])
plot(as.vector(rankPathNetWithoutNull[Endometrial_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 110), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankPathNetWithNull[Endometrial_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Endometrial_cancer]), labels=dataSets[Endometrial_cancer], cex.axis=1)
title("Endometrial cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

dev.off()


############################# BIAS ############################# 
path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
dataSets=c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
           "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")

TargetPW = c("Renal cell carcinoma", "Alzheimers disease", "Thyroid cancer", "Dilated cardiomyopathy",
             "Thyroid cancer", "Colorectal cancer", "Alzheimers disease", "Alzheimers disease",
             "Alzheimers disease", "Prostate cancer", "Prostate cancer", "Colorectal cancer", 
             "Huntingtons disease", "Colorectal cancer", "Acute myeloid leukemia", 
             "Renal cell carcinoma", "Pancreatic cancer", "Pancreatic cancer", "Non-small cell lung cancer",
             "Non-small cell lung cancer", "Glioma", "Parkinsons disease", "Parkinsons disease", "Glioma",
             "Acute myeloid leukemia", "Acute myeloid leukemia", "Alzheimers disease", 
             "Type II diabetes mellitus", "Parkinsons disease", "Colorectal cancer", 
             "Chronic myeloid leukemia", "Chronic myeloid leukemia", "Pancreatic cancer",
             "Colorectal cancer", "Endometrial cancer")

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/PathNet15vs15_2000(2nd run).RData")
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/PathNetTargetPW.RData")


bias0 = vector(mode="numeric", length=nrow(resNullPathNet))
bias1 = vector(mode="numeric", length=nrow(resNullPathNet))
names(bias0) = resNullPathNet$Name
names(bias1) = resNullPathNet$Name

for (i in 1:nrow(resNullPathNet)) {
    # print(i)
    x = resNullPathNet[i,-c(1:3)]
    x =  x[colSums(!is.na(x)) > 0]
    x = as.numeric(x)
    # hist(x, breaks = 50)
    densX = density(x)
    # plot(densX)
    fun = data.frame(x = densX$x, y = densX$y)
    
    # area0 = fun[fun$x <= 0.25 & fun$x > 0 ,]
    # score0 = (area0$y-1)/(area0$x * area0$x)
    # area1 = fun[fun$x >= 0.75 & fun$x < 1 ,]
    # score1 = (area1$y-1)/((1-area1$x) * (1-area1$x))
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

save(dataSets, TargetPW, rankBias0, rankWithNullBias0, rankBias1, rankWithNullBias1, rankWithoutNull, rankWithNull,pValueBias0, 
     pValueWithNullBias0,pValueBias1, pValueWithNullBias1,pvalueWithoutNull, pvalueWithNull, file = "PathNet_bias.RData")

####################### DISEASE-Wise ####################### 

boxplot(rankWithoutNull[Renal_cell_carcinoma], rankWithNull[Renal_cell_carcinoma],
        main="Rank(s) of Renal cell carcinoma (PathNet)", 
        ylab=paste0("Rank (", length(Renal_cell_carcinoma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Alzheimers_disease], rankWithNull[Alzheimers_disease],
        main="Rank(s) of Alzheimer's disease (PathNet)", 
        ylab=paste0("Rank (", length(Alzheimers_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Thyroid_cancer], rankWithNull[Thyroid_cancer],
        main="Rank(s) of Thyroid cancer (PathNet)", 
        ylab=paste0("Rank (", length(Thyroid_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Dilated_cardiomyopathy], rankWithNull[Dilated_cardiomyopathy],
        main="Rank(s) of Dilated cardiomyopathy (PathNet)", 
        ylab=paste0("Rank (", length(Dilated_cardiomyopathy), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Colorectal_cancer], rankWithNull[Colorectal_cancer],
        main="Rank(s) of Colorectal cancer (PathNet)", 
        ylab=paste0("Rank (", length(Colorectal_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Prostate_cancer], rankWithNull[Prostate_cancer],
        main="Rank(s) of Prostate cancer (PathNet)", 
        ylab=paste0("Rank (", length(Prostate_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Huntingtons_disease], rankWithNull[Huntingtons_disease],
        main="Rank(s) of Huntington's disease (PathNet)", 
        ylab=paste0("Rank (", length(Huntingtons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Acute_myeloid_leukemia], rankWithNull[Acute_myeloid_leukemia],
        main="Rank(s) of Acute myeloid leukemia (PathNet)", 
        ylab=paste0("Rank (", length(Acute_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Pancreatic_cancer], rankWithNull[Pancreatic_cancer],
        main="Rank(s) of Pancreatic cancer (PathNet)", 
        ylab=paste0("Rank (", length(Pancreatic_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Non_small_cell_lung_cancer], rankWithNull[Non_small_cell_lung_cancer],
        main="Rank(s) of Non-small cell lung cancer (PathNet)", 
        ylab=paste0("Rank (", length(Non_small_cell_lung_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Glioma], rankWithNull[Glioma],
        main="Rank(s) of Glioma (PathNet)", 
        ylab=paste0("Rank (", length(Glioma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Parkinsons_disease], rankWithNull[Parkinsons_disease],
        main="Rank(s) of Parkinson's disease (PathNet)", 
        ylab=paste0("Rank (", length(Parkinsons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Type_II_diabetes_mellitus], rankWithNull[Type_II_diabetes_mellitus],
        main="Rank(s) of Type II diabetes mellitus (PathNet)", 
        ylab=paste0("Rank (", length(Type_II_diabetes_mellitus), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Chronic_myeloid_leukemia], rankWithNull[Chronic_myeloid_leukemia],
        main="Rank(s) of Chronic myeloid leukemia (PathNet)", 
        ylab=paste0("Rank (", length(Chronic_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Endometrial_cancer], rankWithNull[Endometrial_cancer],
        main="Rank(s) of Endometrial cancer (PathNet)", 
        ylab=paste0("Rank (", length(Endometrial_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))

boxplot(rankWithoutNull[Cancer], rankWithNull[Cancer],
        main="Rank(s) of Cancer (PathNet)", 
        ylab=paste0("Rank (", length(Cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull[Non_cancer], rankWithNull[Non_cancer],
        main="Rank(s) of Non-cancer (PathNet)", 
        ylab=paste0("Rank (", length(Non_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(rankWithoutNull, rankWithNull,
        main="Rank(s) overall (PathNet)", 
        ylab=paste0("Rank (", length(dataSets), " datasets)" ), names=c("Without Danube", "With Danube"))


######### PVALUE ##########

boxplot(pvalueWithoutNull[Renal_cell_carcinoma], pvalueWithNull[Renal_cell_carcinoma],
        main="pValue(s) of Renal cell carcinoma (PathNet)", 
        ylab=paste0("pValue (", length(Renal_cell_carcinoma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Alzheimers_disease], pvalueWithNull[Alzheimers_disease],
        main="pValue(s) of Alzheimer's disease (PathNet)", 
        ylab=paste0("pValue (", length(Alzheimers_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Thyroid_cancer], pvalueWithNull[Thyroid_cancer],
        main="pValue(s) of Thyroid cancer (PathNet)", 
        ylab=paste0("pValue (", length(Thyroid_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Dilated_cardiomyopathy], pvalueWithNull[Dilated_cardiomyopathy],
        main="pValue(s) of Dilated cardiomyopathy (PathNet)", 
        ylab=paste0("pValue (", length(Dilated_cardiomyopathy), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Colorectal_cancer], pvalueWithNull[Colorectal_cancer],
        main="pValue(s) of Colorectal cancer (PathNet)", 
        ylab=paste0("pValue (", length(Colorectal_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Prostate_cancer], pvalueWithNull[Prostate_cancer],
        main="pValue(s) of Prostate cancer (PathNet)", 
        ylab=paste0("pValue (", length(Prostate_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Huntingtons_disease], pvalueWithNull[Huntingtons_disease],
        main="pValue(s) of Huntington's disease (PathNet)", 
        ylab=paste0("pValue (", length(Huntingtons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Acute_myeloid_leukemia], pvalueWithNull[Acute_myeloid_leukemia],
        main="pValue(s) of Acute myeloid leukemia (PathNet)", 
        ylab=paste0("pValue (", length(Acute_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Pancreatic_cancer], pvalueWithNull[Pancreatic_cancer],
        main="pValue(s) of Pancreatic cancer (PathNet)", 
        ylab=paste0("pValue (", length(Pancreatic_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Non_small_cell_lung_cancer], pvalueWithNull[Non_small_cell_lung_cancer],
        main="pValue(s) of Non-small cell lung cancer (PathNet)", 
        ylab=paste0("pValue (", length(Non_small_cell_lung_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Glioma], pvalueWithNull[Glioma],
        main="pValue(s) of Glioma (PathNet)", 
        ylab=paste0("pValue (", length(Glioma), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Parkinsons_disease], pvalueWithNull[Parkinsons_disease],
        main="pValue(s) of Parkinson's disease (PathNet)", 
        ylab=paste0("pValue (", length(Parkinsons_disease), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Type_II_diabetes_mellitus], pvalueWithNull[Type_II_diabetes_mellitus],
        main="pValue(s) of Type II diabetes mellitus (PathNet)", 
        ylab=paste0("pValue (", length(Type_II_diabetes_mellitus), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Chronic_myeloid_leukemia], pvalueWithNull[Chronic_myeloid_leukemia],
        main="pValue(s) of Chronic myeloid leukemia (PathNet)", 
        ylab=paste0("pValue (", length(Chronic_myeloid_leukemia), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Endometrial_cancer], pvalueWithNull[Endometrial_cancer],
        main="pValue(s) of Endometrial cancer (PathNet)", 
        ylab=paste0("pValue (", length(Endometrial_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))

boxplot(pvalueWithoutNull[Cancer], pvalueWithNull[Cancer],
        main="pValue(s) of Cancer (PathNet)", 
        ylab=paste0("pValue (", length(Cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull[Non_cancer], pvalueWithNull[Non_cancer],
        main="pValue(s) of Non-cancer (PathNet)", 
        ylab=paste0("pValue (", length(Non_cancer), " datasets)" ), names=c("Without Danube", "With Danube"))
boxplot(pvalueWithoutNull, pvalueWithNull,
        main="pValue(s) overall (PathNet)", 
        ylab=paste0("pValue (", length(dataSets), " datasets)" ), names=c("Without Danube", "With Danube"))


################################# Pearson's Bias ################################# 
library(sfsmisc)
path = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("PathNet15vs15_2000(2nd run).RData")
load("PathNetTargetPW.RData")
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


bias = vector(mode="numeric", length=nrow(resNullPathNet))
names(bias) = resNullPathNet$Name

for (i in 1:nrow(resNullPathNet)) {
    x = as.numeric(resNullPathNet[i,-c(1:3)])
    x =  na.omit(x)
    bias[i] = Pearson(x)
}

bias0List = bias[which(bias >= 0.1)]
bias1List = bias[which(bias <= -0.1)]

bias0List = names(bias0List)
bias1List = names(bias1List)
bias1List[which(bias1List == "Alzheimers disease")] = "Alzheimer's disease"
bias1List[which(bias1List == "Parkinsons disease")] = "Parkinson's disease"
bias1List[which(bias1List == "Huntingtons disease")] = "Huntington's disease"


save(bias0List, bias1List, file = "PathNet_bias.RData")

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
#      file = "PathNet_bias.RData")
