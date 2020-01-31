#source("http://bioconductor.org/biocLite.R")
#biocLite("CePa")
library("ROntoTools")

#setwd("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/")
#path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/"
setwd("/wsu/home/gd/gd03/gd0393/Pathway/dataset/");
path="/wsu/home/gd/gd03/gd0393/Pathway/dataset/"

iterno=2000
PVAL <- 0.05
maxDE=400

#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/KEGG65.150Pathways.RData")
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
#samePlat=c("GSE5281_EC", "GSE5281_HIP","GSE5281_VCX")

#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/mygeneSym.RData")
#load("/wsu/home/gd/gd03/gd0393/Pathway/mygeneSym.RData")

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
    #print(kpn[pw])
    DEinPW <- length(intersect(keggnodes[[pw]],names(DEGenes)))
    if (DEinPW == 0) { temp[pw,"pvalue"] <- NA } 
    else {
      DEnotInPW = length(DEGenes) - DEinPW
      notDEinPW = length(keggnodes[[pw]]) - DEinPW
      notDEnotinPW = length(pvalues) - length(DEGenes) - length(keggnodes[[pw]]) + DEinPW
      table <- matrix(c(DEinPW, notDEinPW, DEnotInPW, notDEnotinPW), nrow = 2, dimnames = list(Genes = c("DE", "Not DE"),Pathway = c("in PW", "not in PW")))
      fisher <- fisher.test(table,alternative = "greater")
      #fisher$p.value
      temp[pw,"pvalue"] <- fisher$p.value
    }
  }
  
  temp
}

l <- list()
for (i in 1:loopNo) {
  a=sample(colnames(Samples),nc+nd,replace=FALSE)
  l[[i]] <- a
}

#library(multicore)
resNull <- mclapply(X=l,FUN=generateNull,mc.cores=10)
#resNull <- lapply(X=l,FUN=generateNull)



#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/CePaORA15vs15_400_2000 (original).RData")

noLoop <- length(resNull)
pathName <- resNull[[1]][,"Name"]

resNullFE=data.frame(row.names=pathName, Name=pathName,M=rep(NA,length(pathName)),std=rep(NA,length(pathName)))
for (i in 1:noLoop) {
  if (length(resNull[[i]][,"pvalue"])==0) {resNullFE[,paste("p_",i,sep="")]=NA}
  else {,paste("p_",i,sep="")]=resNull[[i]][,"pvalue"]}
}
resNullFE$M=apply(resNullFE[4:(noLoop+3)],1,FUN=mean,na.rm=T)
resNullFE$std=apply(resNullFE[4:(noLoop+3)],1,FUN=sd,na.rm=T)

# v4: adding: if (DEinPW == 0) { temp[pw,"pvalue"] <- NA } 
file = paste0("FE",nc,"vs",nd,"_400_",loopNo,"_v4.RData",sep="") 
save(resNullFE,l,file=file)


load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/FE15vs15_400_2000_v4.RData")

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullFE_v4.pdf"
pathName <- resNullFE[,"Name"]
pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullFE[i,-c(1,2,3)]), breaks = 50, main = resNullFE[i,]$Name)
}
dev.off()

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullFE_color_v4.pdf"

pdf(file=PDFPath)

for (i in seq(length(pathName)))  {
  print(i)
  hist(as.numeric(resNullFE[i,-c(1,2,3)]), breaks = 50, main = resNullFE[i,]$Name, col = "blue")
}
dev.off()
