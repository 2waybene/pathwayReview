setwd("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/");
path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
# setwd("/wsu/home/gd/gd03/gd0393/Pathway/dataset/");
# path="/wsu/home/gd/gd03/gd0393/Pathway/dataset/"

iterno=5000
PVAL <- 0.05
maxDE=400

library("ROntoTools")
library(SPIA)
load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/KEGG65.150Pathways.RData")
# load("/wsu/home/gd/gd03/gd0393/Pathway/KEGG65.150Pathways.RData")
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
samePlat=c("GSE5281_EC", "GSE5281_HIP","GSE5281_VCX")
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

args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
#nc=args[1]
#nd=args[2]
#loopNo=args[3]
nc=15
nd=15
loopNo=2000

generateNull = function (a) {
  controlDat = Samples[,a[1:nc]]
  
  diseaseDat = Samples[,a[(nc+1):(nc+nd)]]
  
  controlMean <- apply(controlDat,MARGIN=1,mean)
  diseaseMean <- apply(diseaseDat,MARGIN=1,mean)
  
  foldChange <- diseaseMean-controlMean
  names(foldChange) <- paste("hsa:",names(foldChange),sep="")
  
  pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
  
  names(pvalues) <- names(foldChange)
  
  DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
  DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
  DEGenes <- foldChange[names(DEGenes)]
  
  names(DEGenes) <- gsub("hsa:", "", names(DEGenes))
  names(foldChange) <- gsub("hsa:", "", names(foldChange))
  
  res=spia(de=DEGenes, all=names(foldChange), organism="hsa",beta=NULL,nB=2000,plots=FALSE, verbose=TRUE,combine="fisher")
  
  res <- res[,-12]
  rownames(res) <- res[,1]
  res[,1] <- NULL
  res[,dim(res)[2]+1] <- res[,1]
  res[,dim(res)[2]] <- paste0("path:hsa",res[,dim(res)[2]])
  res[,1] <- NULL
  colnames(res)[dim(res)[2]] = "pathway"
  res
}

l <- list()
for (i in 1:loopNo) {
  a=sample(colnames(Samples),nc+nd,replace=FALSE)
  l[[i]] <- a
}

#library(multicore)
# resNull <- mclapply(X=l,FUN=generateNull,mc.cores=20)
resNull <- lapply(X=l,FUN=generateNull)

resNullIA=data.frame(row.names=names(kpn), Name=kpn,M=rep(NA,length(kpn)),std=rep(NA,length(kpn)))
for (i in 1:loopNo) {
  if (length(resNull[[i]]$pG)==0) {resNullIA[,paste("pIA",i,sep="")]=NA}
  else {resNullIA[resNull[[i]]$pathway,paste("pIA",i,sep="")]=resNull[[i]]$pG}
}
resNullIA$M=apply(resNullIA[4:(loopNo+3)],1,FUN=mean,na.rm=T)
resNullIA$std=apply(resNullIA[4:(loopNo+3)],1,FUN=sd,na.rm=T)


file = paste("SPIA",nc,"vs",nd,"_400_",loopNo,".RData",sep="")
save(resNullIA,l,file=file)



load("SPIA15vs15_400_2000.RData")
loopNo = 2000

resNullSPIA=data.frame(row.names=kpn, Name=kpn,M=rep(NA,length(kpn)),std=rep(NA,length(kpn)))
for (i in 1:loopNo) {
  if (length(resNull[[i]]) != 1 ) {
    print(i)
    if (length(resNull[[i]]$pG)==0) {resNullSPIA[,paste("pIA",i,sep="")]=NA}
    else {resNullSPIA[rownames(resNull[[i]]),paste("pIA",i,sep="")]=resNull[[i]]$pG}
  }
}
resNullSPIA$M=apply(resNullSPIA[4:nrow(resNullSPIA)],1,FUN=mean,na.rm=T)
resNullSPIA$std=apply(resNullSPIA[4:nrow(resNullSPIA)],1,FUN=sd,na.rm=T)


save(resNullSPIA,resNull ,file= "SPIA_15vs15_2000.RData")

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullSPIA.pdf"
pathName <- resNullSPIA[,"Name"]
pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullSPIA[i,-c(1,2,3)]), breaks = 0:50/50, main = resNullSPIA[i,]$Name)
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullSPIA_color.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullSPIA[i,-c(1,2,3)]), breaks = 0:50/50, main = resNullSPIA[i,]$Name, col = "blue")
}
dev.off()


PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullSPIA_20bins.pdf"
pathName <- resNullSPIA[,"Name"]
pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullSPIA[i,-c(1,2,3)]), breaks = 0:20/20, main = resNullSPIA[i,]$Name)
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullSPIA_20bins_color.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullSPIA[i,-c(1,2,3)]), breaks = 0:20/20, main = resNullSPIA[i,]$Name, col = "blue")
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullSPIA_density.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  if (!is.na(resNullSPIA[i,3]) ) {
    print(i)
    x <- density(na.omit(as.numeric(resNullSPIA[i,-c(1,2,3)])))
    plot(x, main = resNullSPIA[i,]$Name) 
    polygon(x, col="blue", border="blue")
    # hist(as.numeric(resNullIA[i,-c(1,2,3)]), breaks = 20, main = resNullIA[i,]$Name, col = "blue")
  }
}
dev.off()

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

# load("CePaGSATargetPW.RData")

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/SPIA_15vs15_2000.RData")

pw = pwNames$Name
resNullSPIA = resNullSPIA[intersect(pw,rownames(resNullSPIA)),]

bias = vector(mode="numeric", length=nrow(resNullSPIA))
for(i in 1:nrow(resNullSPIA)) {
  bias[i] = NA
}
names(bias) = resNullSPIA$Name

for (i in 1:nrow(resNullSPIA)) {
  if (!is.na(resNullSPIA[i,3]) ) {
    print(i)
    x = as.numeric(resNullSPIA[i,-c(1:3)])
    x =  na.omit(x)
    bias[i] = Pearson(x)
  }
}

bias0List = bias[which(bias >= 0.1)]
bias1List = bias[which(bias <= -0.1)]

bias0List = names(bias0List)
bias1List = names(bias1List)
bias0List %in% pwNames$Name
bias1List %in% pwNames$Name

save(bias0List, bias1List, file = "SPIA_bias.RData")

