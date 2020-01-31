#setwd("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/");
#path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/"
setwd("/wsu/home/gd/gd03/gd0393/Pathway/dataset/");
path="/wsu/home/gd/gd03/gd0393/Pathway/dataset/"
#dataSets=c("GSE1297","GSE28146","GSE5281_EC","GSE5281_HIP","GSE5281_MTG","GSE5281_PC", "GSE5281_SFG", "GSE5281_VCX")
dataSets=c("GSE6956_AA", "GSE6956_C")
iterno=1000
targetGeneSets="05010"

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/PADOG_plus2_15vs15_400_final.RData")
PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullPADOG.pdf"
pathName <- resNullAll[,"Name"]
pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullAll[i,-c(1,2,3)]), breaks = 50, main = resNullAll[i,]$Name)
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullPADOG_color.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullAll[i,-c(1,2,3)]), breaks = 50, main = resNullAll[i,]$Name, col = "blue")
}
dev.off()


load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/PADOG_plus2_15vs15_400_final.RData")
PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullPADOG_20bins.pdf"
pathName <- resNullAll[,"Name"]
pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullAll[i,-c(1,2,3)]), breaks = 20, main = resNullAll[i,]$Name)
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullPADOG_20bins_color.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullAll[i,-c(1,2,3)]), breaks = 20, main = resNullAll[i,]$Name, col = "blue")
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullPADOG_density.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
    print(i)
    x <- density(na.omit(as.numeric(resNullAll[i,-c(1,2,3)])))
    plot(x, main = resNullAll[i,]$Name) 
    polygon(x, col="blue", border="blue")
    # hist(as.numeric(resNullIA[i,-c(1,2,3)]), breaks = 20, main = resNullIA[i,]$Name, col = "blue")
  
}
dev.off()

