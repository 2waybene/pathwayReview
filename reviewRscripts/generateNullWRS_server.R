#source("http://bioconductor.org/biocLite.R")
#biocLite("CePa")
library("ROntoTools")

# setwd("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/")
# path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/"
setwd("/wsu/home/gd/gd03/gd0393/Pathway/dataset/");
path="/wsu/home/gd/gd03/gd0393/Pathway/dataset/"

iterno=2000
PVAL <- 0.05
maxDE=400

# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/KEGG65.150Pathways.RData")
load("/wsu/home/gd/gd03/gd0393/Pathway/KEGG65.150Pathways.RData")
kpg <- setEdgeWeights(kpg)
kpg <- setNodeWeights(kpg, defaultWeight = 1)

getnodes=function(x){
  return (x@nodes);}


#generate the null distribution
samePlat=c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
           "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")
# samePlat=c("GSE5281_EC", "GSE5281_HIP","GSE5281_VCX")

# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/mygeneSym.RData")
# load("/wsu/home/gd/gd03/gd0393/Pathway/mygeneSym.RData")

Samples <- NULL
intersectGenes <- NULL
for (i in 1:length(samePlat)) {
  dataset=samePlat[i]
  print(dataset)
  filename=paste(path,dataset,"/",dataset,"_annotation.txt",sep="")
  annotation=as.character((read.table(filename,header=F,sep="\t",stringsAsFactors=F))[1,1])
  keggnodes=lapply(kpg,getnodes)
  allGenes=unique(unlist(keggnodes))
  ano=paste(substring(annotation, 1, nchar(annotation)-3),"ENTREZID",sep="")
  require(annotation, character.only = TRUE)
  geneAno=as.data.frame(get(ano))
  geneAnno=geneAno[,2]
  names(geneAnno)=geneAno[,1]
  allProbes=unique(geneAno$probe_id)
  
  load(paste(path,dataset,"/",dataset,".RData",sep=""))      
  data <- get(paste("data_",dataset,sep=""))
  
  group <- get(paste("group_",dataset,sep=""))
  cn = rownames(group[group$Group=="c",])
  data <- data[,cn]
  
  data <- 2^data
  data <- data[rownames(data)%in%allProbes,]
  rownames(data) <- as.character(geneAnno[rownames(data)])
  data <- aggregate(data,by=list(rownames(data)),FUN=median)
  #data <- data[paste("hsa:",data$Group.1,sep="")%in%allGenes,]
  data <- data[paste("hsa:",data$Group.1,sep="")%in%allGenes,]
  rownames(data)=data$Group.1
  data <- data[,!colnames(data)%in%c("Group.1")]
  
  
  
  if (is.null(intersectGenes)) {
    intersectGenes = rownames(data)
  } else {
    intersectGenes = intersect(intersectGenes, rownames(data))
  }
  if (is.null(Samples)) {
    Samples = data
  } else {
    Samples <- cbind(Samples[intersectGenes,],data[intersectGenes,])
  }
  
}

#args <- commandArgs(trailingOnly = TRUE)
#args <- as.numeric(args)
#nc=args[1]
#nd=args[2]
#loopNo=args[3]
nc=15
nd=15
#loopNo=5000
loopNo=2000

generateNull = function (a) {
  controlDat = Samples[,a[1:nc]]
  
  diseaseDat = Samples[,a[(nc+1):(nc+nd)]]
  #wholeDat <- cbind(controlDat, diseaseDat)
  ##########
  
  controlMean <- log(apply(controlDat,MARGIN=1,mean),base=2)
  diseaseMean <- log(apply(diseaseDat,MARGIN=1,mean),base=2)
  
  foldChange <- diseaseMean-controlMean
  names(foldChange) <- paste("hsa:",names(foldChange),sep="")
  
  pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
  
  names(pvalues) <- names(foldChange)
  
  DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
  DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
  DEGenes <- foldChange[names(DEGenes)]
  
  temp=data.frame(row.names=names(kpn), Name=kpn,pvalue=rep(NA,length(kpn)))
  for (pw in 1:length(kpn)) {
    DEhit <- DEGenes[intersect(keggnodes[[pw]],names(DEGenes))]
    DEmiss <- DEGenes[!DEGenes%in%DEhit]
    
    
    if (length(DEhit) == 0 || length(DEmiss) == 0) { pval = 1;}
    else {
      wrs <- wilcox.test(DEhit,DEmiss)
      pval = wrs$p.value
    }
    
    temp[names(kpn)[pw],"pvalue"] <- pval
  }
  
  temp
}

l <- list()
for (i in 1:loopNo) {
  a=sample(colnames(Samples),nc+nd,replace=FALSE)
  l[[i]] <- a
}

#library(multicore)
resNull <- mclapply(X=l,FUN=generateNull,mc.cores=20)
# resNull <- lapply(X=l,FUN=generateNull)

#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/CePaORA15vs15_400_2000 (original).RData")

noLoop <- length(resNull)
pathName <- resNull[[1]][,"Name"]

resNullWRS=data.frame(row.names=pathName, Name=pathName,M=rep(NA,length(pathName)),std=rep(NA,length(pathName)))
for (i in 1:noLoop) {
  if (length(resNull[[i]][,"pvalue"])==0) {resNullWRS[,paste("p_",i,sep="")]=NA}
  else {resNullWRS[as.vector(resNull[[i]][,"Name"]),paste("p_",i,sep="")]=resNull[[i]][,"pvalue"]}
}
resNullWRS$M=apply(resNullWRS[4:(noLoop+3)],1,FUN=mean,na.rm=T)
resNullWRS$std=apply(resNullWRS[4:(noLoop+3)],1,FUN=sd,na.rm=T)

file = paste0("WRS",nc,"vs",nd,"_400_",loopNo,".RData",sep="") 
save(resNullWRS,l,file=file)


load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/WRS15vs15_400_2000(2nd run).RData")
PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullWRS (2nd run).pdf"
pathName <- resNullWRS[,"Name"]
pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullWRS[i,-c(1,2,3)]), breaks = 50, main = resNullWRS[i,]$Name)
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullWRS_color (2nd run).pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullWRS[i,-c(1,2,3)]), breaks = 50, main = resNullWRS[i,]$Name, col = "blue")
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullWRS_20bins.pdf"
pathName <- resNullWRS[,"Name"]
pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullWRS[i,-c(1,2,3)]), breaks = 20, main = resNullWRS[i,]$Name)
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullWRS_20bins_color.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullWRS[i,-c(1,2,3)]), breaks = 20, main = resNullWRS[i,]$Name, col = "blue")
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullWRS_density.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  x <- density(na.omit(as.numeric(resNullWRS[i,-c(1,2,3)])))
  plot(x, main = resNullWRS[i,]$Name) 
  polygon(x, col="blue", border="blue")
  # hist(as.numeric(resNullIA[i,-c(1,2,3)]), breaks = 20, main = resNullIA[i,]$Name, col = "blue")
  
}
dev.off()
