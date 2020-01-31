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

rankWithoutNull <- rep(NA, 35)
pvalueWithoutNull <- rep(NA, 35)

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsGSEA.RData")

for(i in 1:length(Renal_cell_carcinoma)) {
  temp = as.data.frame(myData[[Renal_cell_carcinoma[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Renal_cell_carcinoma[i]] <- temp["path:hsa05211","rank"]
  pvalueWithoutNull[Renal_cell_carcinoma[i]] <- p["path:hsa05211"]
}

for(i in 1:length(Alzheimers_disease)) {
  temp = as.data.frame(myData[[Alzheimers_disease[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Alzheimers_disease[i]] <- temp["path:hsa05010","rank"]
  pvalueWithoutNull[Alzheimers_disease[i]] <- p["path:hsa05010"]
}

for(i in 1:length(Thyroid_cancer)) {
  temp = as.data.frame(myData[[Thyroid_cancer[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Thyroid_cancer[i]] <- temp["path:hsa05216","rank"]
  pvalueWithoutNull[Thyroid_cancer[i]] <- p["path:hsa05216"]
}

for(i in 1:length(Colorectal_cancer)) {
  temp = as.data.frame(myData[[Colorectal_cancer[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Colorectal_cancer[i]] <- temp["path:hsa05210","rank"]
  pvalueWithoutNull[Colorectal_cancer[i]] <- p["path:hsa05210"]
}

for(i in 1:length(Prostate_cancer)) {
  temp = as.data.frame(myData[[Prostate_cancer[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Prostate_cancer[i]] <- temp["path:hsa05215","rank"]
  pvalueWithoutNull[Prostate_cancer[i]] <- p["path:hsa05215"]
}

for(i in 1:length(Acute_myeloid_leukemia)) {
  temp = as.data.frame(myData[[Acute_myeloid_leukemia[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Acute_myeloid_leukemia[i]] <- temp["path:hsa05221","rank"]
  pvalueWithoutNull[Acute_myeloid_leukemia[i]] <- p["path:hsa05221"]
}

for(i in 1:length(Pancreatic_cancer)) {
  temp = as.data.frame(myData[[Pancreatic_cancer[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Pancreatic_cancer[i]] <- temp["path:hsa05212","rank"]
  pvalueWithoutNull[Pancreatic_cancer[i]] <- p["path:hsa05212"]
}

for(i in 1:length(Non_small_cell_lung_cancer)) {
  temp = as.data.frame(myData[[Non_small_cell_lung_cancer[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Non_small_cell_lung_cancer[i]] <- temp["path:hsa05223","rank"]
  pvalueWithoutNull[Non_small_cell_lung_cancer[i]] <- p["path:hsa05223"]
}

for(i in 1:length(Glioma)) {
  temp = as.data.frame(myData[[Glioma[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Glioma[i]] <- temp["path:hsa05214","rank"]
  pvalueWithoutNull[Glioma[i]] <- p["path:hsa05214"]
}

for(i in 1:length(Parkinsons_disease)) {
  temp = as.data.frame(myData[[Parkinsons_disease[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Parkinsons_disease[i]] <- temp["path:hsa05012","rank"]
  pvalueWithoutNull[Parkinsons_disease[i]] <- p["path:hsa05012"]
}

for(i in 1:length(Chronic_myeloid_leukemia)) {
  temp = as.data.frame(myData[[Chronic_myeloid_leukemia[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Chronic_myeloid_leukemia[i]] <- temp["path:hsa05012","rank"]
  pvalueWithoutNull[Chronic_myeloid_leukemia[i]] <- p["path:hsa05012"]  
}

for(i in 1:length(Type_II_diabetes_mellitus)) {
  temp = as.data.frame(myData[[Type_II_diabetes_mellitus[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Type_II_diabetes_mellitus[i]] <- temp["path:hsa04930","rank"]
  pvalueWithoutNull[Type_II_diabetes_mellitus[i]] <- p["path:hsa04930"] 
}

for(i in 1:length(Dilated_cardiomyopathy)) {
  temp = as.data.frame(myData[[Dilated_cardiomyopathy[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Dilated_cardiomyopathy[i]] <- temp["path:hsa05414","rank"]
  pvalueWithoutNull[Dilated_cardiomyopathy[i]] <- p["path:hsa05414"] 
}

for(i in 1:length(Huntingtons_disease)) {
  temp = as.data.frame(myData[[Huntingtons_disease[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Huntingtons_disease[i]] <- temp["path:hsa05016","rank"]
  pvalueWithoutNull[Huntingtons_disease[i]] <- p["path:hsa05016"] 
}

for(i in 1:length(Endometrial_cancer)) {
  temp = as.data.frame(myData[[Endometrial_cancer[i]]])
  RowNames <- rownames(temp)
  p <- as.numeric(levels(temp[,"pvalue"]))[temp[,"pvalue"]]
  names(p) <- RowNames
  rankWithoutNull[Endometrial_cancer[i]] <- temp["path:hsa05213","rank"]
  pvalueWithoutNull[Endometrial_cancer[i]] <- p["path:hsa05213"] 
}

save(rankWithoutNull, pvalueWithoutNull, file=paste0(path,"GSEATargetPW.RData"))


boxplot(as.vector(as.numeric(rankWithoutNull)),main ="Ranks of target pathways", ylim=c(0, 120))
boxplot(as.vector(as.numeric(pvalueWithoutNull)),main ="p-values of target pathways", ylim = c(0,1))


vioplot(as.vector(as.numeric(rankWithoutNull)), col = "grey", ylim = c(0,120))
title("GSEA without DANUBE", ylab = "Rank of target pathways")
median(as.vector(as.numeric(rankWithoutNull)))

vioplot(as.vector(as.numeric(pvalueWithoutNull)), col = "grey", ylim = c(0,1))
title("GSEA without DANUBE", ylab = "pValie of target pathways")
median(as.vector(as.numeric(pvalueWithoutNull)))


vioplot(as.vector(as.numeric(rankWSRWithNull)), col = "grey", ylim = c(0,120))
title("WRS with DANUBE", ylab = "Rank of target pathways")
median(as.numeric(rankWSRWithNull))
vioplot(as.vector(as.numeric(pValueWSRWithNull)), col = "grey", ylim = c(0,1))
title("WRS with DANUBE", ylab = "pValue of target pathways")
median(as.vector(as.numeric(pValueWSRWithNull)))


col <- c("orange", "blue")

as.vector(rankWithoutNull[Renal_cell_carcinoma])
as.vector(rankWSRWithNull[Renal_cell_carcinoma])
plot(as.vector(rankWithoutNull[Renal_cell_carcinoma]), col = col[1], type = 'b', 
     ylim=c(0, 130), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Renal_cell_carcinoma]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Renal_cell_carcinoma]), labels=dataSets[Renal_cell_carcinoma], cex.axis=1)
title("Renal cell carcinoma")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Alzheimers_disease])
as.vector(rankWSRWithNull[Alzheimers_disease])
plot(as.vector(rankWithoutNull[Alzheimers_disease]), col = col[1], type = 'b', 
     ylim=c(0, 120), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Alzheimers_disease]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Alzheimers_disease]), labels=dataSets[Alzheimers_disease], cex.axis=1)
title("Alzheimer's disease")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Thyroid_cancer])
as.vector(rankWSRWithNull[Thyroid_cancer])
plot(as.vector(rankWithoutNull[Thyroid_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 110), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Thyroid_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Thyroid_cancer]), labels=dataSets[Thyroid_cancer], cex.axis=1)
title("Thyroid cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Dilated_cardiomyopathy])
as.vector(rankWSRWithNull[Dilated_cardiomyopathy])
plot(as.vector(rankWithoutNull[Dilated_cardiomyopathy]), col = col[1], type = 'b', 
     ylim=c(0, 90), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Dilated_cardiomyopathy]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Dilated_cardiomyopathy]), labels=dataSets[Dilated_cardiomyopathy], cex.axis=1)
title("Dilated cardiomyopathy")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Colorectal_cancer])
as.vector(rankWSRWithNull[Colorectal_cancer])
plot(as.vector(rankWithoutNull[Colorectal_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 90), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Colorectal_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Colorectal_cancer]), labels=dataSets[Colorectal_cancer], cex.axis=1)
title("Colorectal cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Prostate_cancer])
as.vector(rankWSRWithNull[Prostate_cancer])
plot(as.vector(rankWithoutNull[Prostate_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 80), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Prostate_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Prostate_cancer]), labels=dataSets[Prostate_cancer], cex.axis=1)
title("Prostate cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Huntingtons_disease])
as.vector(rankWSRWithNull[Huntingtons_disease])
plot(as.vector(rankWithoutNull[Huntingtons_disease]), col = col[1], type = 'b', 
     ylim=c(0, 20), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Huntingtons_disease]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Huntingtons_disease]), labels=dataSets[Huntingtons_disease], cex.axis=1)
title("Huntington's disease")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Acute_myeloid_leukemia])
as.vector(rankWSRWithNull[Acute_myeloid_leukemia])
plot(as.vector(rankWithoutNull[Acute_myeloid_leukemia]), col = col[1], type = 'b', 
     ylim=c(0, 50), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Acute_myeloid_leukemia]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Acute_myeloid_leukemia]), labels=dataSets[Acute_myeloid_leukemia], cex.axis=1)
title("Acute myeloid leukemia")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Pancreatic_cancer])
as.vector(rankWSRWithNull[Pancreatic_cancer])
plot(as.vector(rankWithoutNull[Pancreatic_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 25), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Pancreatic_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Pancreatic_cancer]), labels=dataSets[Pancreatic_cancer], cex.axis=1)
title("Pancreatic cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Non_small_cell_lung_cancer])
as.vector(rankWSRWithNull[Non_small_cell_lung_cancer])
plot(as.vector(rankWithoutNull[Non_small_cell_lung_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 50), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Non_small_cell_lung_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Non_small_cell_lung_cancer]), labels=dataSets[Non_small_cell_lung_cancer], cex.axis=1)
title("Non-small cell lung cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Glioma])
as.vector(rankWSRWithNull[Glioma])
plot(as.vector(rankWithoutNull[Glioma]), col = col[1], type = 'b', 
     ylim=c(0, 20), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Glioma]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Glioma]), labels=dataSets[Glioma], cex.axis=1)
title("Glioma")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Parkinsons_disease])
as.vector(rankWSRWithNull[Parkinsons_disease])
plot(as.vector(rankWithoutNull[Parkinsons_disease]), col = col[1], type = 'b', 
     ylim=c(0, 100), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Parkinsons_disease]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Parkinsons_disease]), labels=dataSets[Parkinsons_disease], cex.axis=1)
title("Parkinson's disease")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Type_II_diabetes_mellitus])
as.vector(rankWSRWithNull[Type_II_diabetes_mellitus])
plot(as.vector(rankWithoutNull[Type_II_diabetes_mellitus]), col = col[1], type = 'b', 
     ylim=c(0, 25), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Type_II_diabetes_mellitus]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Type_II_diabetes_mellitus]), labels=dataSets[Type_II_diabetes_mellitus], cex.axis=1)
title("Type II diabetes mellitus")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Chronic_myeloid_leukemia])
as.vector(rankWSRWithNull[Chronic_myeloid_leukemia])
plot(as.vector(rankWithoutNull[Chronic_myeloid_leukemia]), col = col[1], type = 'b', 
     ylim=c(0, 95), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Chronic_myeloid_leukemia]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Chronic_myeloid_leukemia]), labels=dataSets[Chronic_myeloid_leukemia], cex.axis=1)
title("Chronic myeloid leukemia")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

as.vector(rankWithoutNull[Endometrial_cancer])
as.vector(rankWSRWithNull[Endometrial_cancer])
plot(as.vector(rankWithoutNull[Endometrial_cancer]), col = col[1], type = 'b', 
     ylim=c(0, 110), xaxt='n', las=2, xlab = '', cex.axis=1, ylab = 'Target pathway rank', pch = 0)
lines(as.vector(rankWSRWithNull[Endometrial_cancer]), col = col[2], type = 'b')
axis(1, at=1:length(dataSets[Endometrial_cancer]), labels=dataSets[Endometrial_cancer], cex.axis=1)
title("Endometrial cancer")
legend('topright', col=col, lwd=c(1, 1.5, 2), legend=c('without adjustment','with adjustment'), pch = c(0,1))

dev.off()


######################## BIAS ####################################

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


load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/WRS15vs15_400_2000(2nd run).RData")
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/WSRTargetPathway.RData")
rankWithoutNull = as.numeric(rankWithoutNull)
pvalueWithoutNull = as.numeric(pvalueWithoutNull)
dim(resNullWRS)

bias0 = vector(mode="numeric", length=nrow(resNullWRS))
bias1 = vector(mode="numeric", length=nrow(resNullWRS))
names(bias0) = resNullWRS$Name
names(bias1) = resNullWRS$Name

for (i in 1:nrow(resNullWRS)) {
  
  # print(i)
  x = resNullWRS[i,-c(1:3)]
  x =  x[colSums(!is.na(x)) > 0]
  x = as.numeric(x)
  # hist(x, breaks = 50)
  densX = density(x)
  # plot(densX)
  fun = data.frame(x = densX$x, y = densX$y)
  bias0[i] = max(fun[fun$x <= 0.5,2] - 1)
  bias1[i] = max(fun[fun$x > 0.5,2] - 1)
  
}


# topBias = c(quantile(bias0)[4], quantile(bias1)[4])
# bias0pw = names(bias0[which(bias0 >= topBias[1])])
# bias1pw = names(bias1[which(bias1 >= topBias[2])])

bias0List = bias0[TargetPW]
bias1List = bias1[TargetPW]
topBias = c(quantile(bias0List)[4], quantile(bias1List)[4])
bias0List =which(bias0List >= topBias[1])
bias1List =which(bias1List >= topBias[2])

rankBias0 = rankWithoutNull[bias0List]
rankBias1 = rankWithoutNull[bias1List]
pValueBias0 = pvalueWithoutNull[bias0List]
pValueBias1 = pvalueWithoutNull[bias1List]

vioplot(as.vector(rankBias0), col = "grey", ylim=c(0, 120))
title("Ranks of target pathways (biased towards 0)")
median(rankBias0)
vioplot(as.vector(pValueBias0), col = "grey",ylim=c(0, 1))
title("pValue of target pathways (biased towards 0)")
median(pValueBias0)


vioplot(as.vector(rankBias1), col = "grey",ylim=c(0, 120))
title("Ranks of target pathways (biased towards 1)")
median(rankBias1)
vioplot(as.vector(pValueBias1), col = "grey", ylim=c(0, 1))
title("pValue of target pathways (biased towards 1)")
median(pValueBias1)

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/FEWithNull.RData")

rankWithNullBias0 = rankWSRWithNull[bias0List]
rankWithNullBias1 = rankWSRWithNull[bias1List]
pValueWithNullBias0 = pValueWSRWithNull[bias0List]
pValueWithNullBias1 = pValueWSRWithNull[bias1List]



vioplot(as.vector(rankWithNullBias0), col = "grey", ylim=c(0, 120))
title("Ranks of target pathways (biased towards 0)")
median(rankWithNullBias0)
vioplot(as.vector(pValueWithNullBias0), col = "grey",ylim=c(0, 1))
title("pValue of target pathways (biased towards 0)")
median(pValueWithNullBias0)


vioplot(as.vector(rankWithNullBias1), col = "grey",ylim=c(0, 120))
title("Ranks of target pathways (biased towards 1)")
median(rankWithNullBias1)
vioplot(as.vector(pValueWithNullBias1), col = "grey", ylim=c(0, 1))
title("pValue of target pathways (biased towards 1)")
median(pValueWithNullBias1)

###################################### Pearson's Bias ################################# 

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

load("GSEATargetPW.RData")

load("GSEA15vs15_400_2000.RData")
# save(dataSets, TargetPW, rankWithNull, rankWithoutNull, pvalueWithNull, pvalueWithoutNull, file = "ROntoToolsTargetPW.RData")

# dim(resNullGSEA)

bias = vector(mode="numeric", length=nrow(resNullGSEA))
names(bias) = resNullGSEA$Name

for (i in 1:nrow(resNullGSEA)) {
  
  x = unlist(resNullGSEA[i,-c(1,2,3)]) 
  x = as.numeric(as.character(temp))
  x =  na.omit(x)
  bias[i] = Pearson(x)
}

bias0List = bias[which(bias >= 0.1)]
bias1List = bias[which(bias <= -0.1)]

bias0List = names(bias0List)
bias1List = names(bias1List)
bias0List %in% pwNames$Name
bias1List %in% pwNames$Name

save(bias0List, bias1List, file = "GSEA_bias.RData")

