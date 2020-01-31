#setwd("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/");
#path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/"
setwd("/wsu/home/gd/gd03/gd0393/Pathway/dataset/");
path="/wsu/home/gd/gd03/gd0393/Pathway/dataset/"
#dataSets=c("GSE1297","GSE28146","GSE5281_EC","GSE5281_HIP","GSE5281_MTG","GSE5281_PC", "GSE5281_SFG", "GSE5281_VCX")
dataSets=c("GSE3467","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE8671","GSE8762",
           "GSE9348","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")
iterno=1000


load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/GSA_plus2_15vs15_400_final.RData")
PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullGSA.pdf"
pathName <- resNullAll[,"Name"]
pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullAll[i,-c(1,2,3)]), breaks = 50, main = resNullAll[i,]$Name)
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullGSA_color.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullAll[i,-c(1,2,3)]), breaks = 50, main = resNullAll[i,]$Name, col = "blue")
}
dev.off()


load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/GSA_plus2_15vs15_400_final.RData")
PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullGSA_20bins.pdf"
pathName <- resNullAll[,"Name"]
pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullAll[i,-c(1,2,3)]), breaks = 20, main = resNullAll[i,]$Name)
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullGSA_20bins_color.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullAll[i,-c(1,2,3)]), breaks = 20, main = resNullAll[i,]$Name, col = "blue")
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullGSA_density.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
    print(i)
    x <- density(na.omit(as.numeric(resNullAll[i,-c(1,2,3)])))
    plot(x, main = resNullAll[i,]$Name) 
    polygon(x, col="blue", border="blue")
    # hist(as.numeric(resNullIA[i,-c(1,2,3)]), breaks = 20, main = resNullIA[i,]$Name, col = "blue")
  
}
dev.off()

