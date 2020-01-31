#source("http://bioconductor.org/biocLite.R")
#biocLite("CePa")
# library("ROntoTools")
library(CePa)
data(PID.db)

library(org.Hs.eg.db)

# setwd("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/");
# path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/"
setwd("/wsu/home/gd/gd03/gd0393/Pathway/dataset/");
path="/wsu/home/gd/gd03/gd0393/Pathway/dataset/"

iterno=1000
PVAL <- 0.05


# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/KEGG65.150Pathways.RData")
#load("/wsu/home/gd/gd03/gd0393/Pathway/KEGG65.150Pathways.RData")
# kpg <- setEdgeWeights(kpg)
# kpg <- setNodeWeights(kpg, defaultWeight = 1)
# 
# getnodes=function(x){
#   return (x@nodes);}


#generate the null distribution
samePlat=c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
           "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")
# samePlat=c("GSE781", "GSE1297","GSE3467","GSE3585")

# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/mygeneSym.RData")
load("/wsu/home/gd/gd03/gd0393/Pathway/mygeneSym.RData")

Samples <- NULL
intersectGenes <- NULL
for (i in 1:length(samePlat)) {
  dataset=samePlat[i]
  print(dataset)
  filename=paste(path,dataset,"/",dataset,"_annotation.txt",sep="")
  annotation=as.character((read.table(filename,header=F,sep="\t",stringsAsFactors=F))[1,1])
  # keggnodes=lapply(kpg,getnodes)
  # allGenes=unique(unlist(keggnodes))
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
  # data <- data[paste("hsa:",data$Group.1,sep="")%in%allGenes,]
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

maxDE=1200
nc=15
nd=15
#loopNo=5000
loopNo=3

generateNull = function (a) {
  controlDat = Samples[,a[1:nc]]
  
  diseaseDat = Samples[,a[(nc+1):(nc+nd)]]
  
  controlMean <- apply(controlDat,MARGIN=1,mean)
  diseaseMean <- apply(diseaseDat,MARGIN=1,mean)
  
  foldChange <- diseaseMean-controlMean
  names(foldChange) <- paste("hsa:",names(foldChange),sep="")
  
  
  pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
  
  names(pvalues) <- names(foldChange)
  
  allgeneID=gsub(pattern = "hsa:", replacement = "", x = names(foldChange))
  allgeneSYM=as.character(mygeneSym[allgeneID])
  
  DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
  DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
  DEGenes <- foldChange[names(DEGenes)]
  
  myDEgeneIDS=gsub(pattern = "hsa:", replacement = "", x = names(DEGenes))
  myDEgeneSYM=as.character(mygeneSym[myDEgeneIDS])
  
  length(myDEgeneSYM)
  length(allgeneSYM)
  
  res.ora1 =cepa.all(dif=myDEgeneSYM, bk=allgeneSYM, pc=PID.db$KEGG)
  temp <- p.table(res.ora1)
  min <- apply(temp, 1, FUN=min)
  temp <- data.frame(Pathway = row.names(temp), pValue = min)
  
  temp
}

l <- list()
for (i in 1:loopNo) {
  a=sample(colnames(Samples),nc+nd,replace=FALSE)
  l[[i]] <- a
}

resNull <- mclapply(X=l,FUN=generateNull,mc.cores=10)
# resNull <- lapply(X=l,FUN=generateNull)

noLoop <- length(resNull)
pathName <- resNull[[1]][,"Pathway"]

resNullCePaORA=data.frame(row.names=pathName, Name=pathName,M=rep(NA,length(pathName)),std=rep(NA,length(pathName)))
for (i in 1:noLoop) {
  # if (length(resNull[[i]][,"pValue"])==0) {
  #   resNullCePaORA[,paste("p_",i,sep="")]=NA
  # } else {
    resNullCePaORA[as.vector(resNull[[i]][,"Pathway"]),paste("p_",i,sep="")]=as.numeric(resNull[[i]][,"pValue"])
    # }
}
resNullCePaORA$M=apply(resNullCePaORA[4:(noLoop+3)],1,FUN=mean,na.rm=T)
resNullCePaORA$std=apply(resNullCePaORA[4:(noLoop+3)],1,FUN=sd,na.rm=T)



library(KEGG.db)
keggid <- "hsa00232"
keggid2keggname <- as.list(KEGGPATHID2NAME)
names(keggid2keggname) <- paste0("hsa",names(keggid2keggname))


x <- as.matrix(keggid2keggname)
x <- matrix(c(row.names(x),x), ncol = 2)
rownames(x) <- x[,1]
colnames(x) <- c("PathwayID", "Pathway")
y <- rownames(resNullCePaORA)
rownames(resNullCePaORA) <- x[rownames(resNullCePaORA),"Pathway"]
resNullCePaORA[,1] <- rownames(resNullCePaORA)
rownames(resNullCePaORA) <- y


file = paste0("CePaORA",nc,"vs",nd,"_400_",loopNo,"_v2.RData",sep="") 
save(resNullCePaORA,file=file)

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/CePaORA15vs15_400_2000_v2.RData")

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullCePaORA.pdf"
pathName <- resNullCePaORA[,"Name"]
pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullCePaORA[i,-c(1,2,3)]), breaks = 0:50/50, main = resNullCePaORA[i,]$Name)
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullCePaORA_color.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullCePaORA[i,-c(1,2,3)]), breaks = 0:50/50, main = resNullCePaORA[i,]$Name, col = "blue")
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullCePaORA_20bins.pdf"
pathName <- resNullCePaORA[,"Name"]
pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullCePaORA[i,-c(1,2,3)]), breaks = 0:20/20, main = resNullCePaORA[i,]$Name)
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullCePaORA_20bins_color.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullCePaORA[i,-c(1,2,3)]), breaks = 0:20/20, main = resNullCePaORA[i,]$Name, col = "blue")
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullCePaORA_density.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  x <- density(na.omit(as.numeric(resNullCePaORA[i,-c(1,2,3)])))
  plot(x, main = resNullCePaORA[i,]$Name) 
  polygon(x, col="blue", border="blue")
  # hist(as.numeric(resNullIA[i,-c(1,2,3)]), breaks = 20, main = resNullIA[i,]$Name, col = "blue")
  
}
dev.off()

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/ROntoTools15vs15_400_5000_onlyHealthy.RData")
pw = rownames(resNullIA)
pw = gsub("path:","", pw)
resNullCePaORA = resNullCePaORA[intersect(pw,rownames(resNullCePaORA)),]
pathName <- resNullCePaORA[,"Name"]
PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/UnderTheNull/generateNullCePaORA_density (intersectKegg150).pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  x <- density(na.omit(as.numeric(resNullCePaORA[i,-c(1,2,3)])))
  plot(x, main = resNullCePaORA[i,]$Name) 
  polygon(x, col="blue", border="blue")
  # hist(as.numeric(resNullIA[i,-c(1,2,3)]), breaks = 20, main = resNullIA[i,]$Name, col = "blue")
  
}
dev.off()


