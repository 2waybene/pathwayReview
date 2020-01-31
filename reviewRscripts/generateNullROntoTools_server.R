#setwd("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/");
#path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/"
setwd("/wsu/home/gd/gd03/gd0393/Pathway/dataset/");
path="/wsu/home/gd/gd03/gd0393/Pathway/dataset/"

iterno=2000
PVAL <- 0.05
maxDE=400

library("ROntoTools")
# kpg <- keggPathwayGraphs("hsa")
# kpg <- setEdgeWeights(kpg)
# kpg <- setNodeWeights(kpg, defaultWeight = 1)
# kpn <- keggPathwayNames("hsa")
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

args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
#nc=args[1]
#nd=args[2]
#loopNo=args[3]
nc=15
nd=15
loopNo=5000

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
  ##########
  #groups = as.vector(c(rep(0, nc), rep(1, nd)))
  #names <- c(colnames(controlDat), colnames(diseaseDat))
  #names(groups) <- names
  #d <- list(cbind(controlDat, diseaseDat), groups)
  
  #resultFC <- logfc(d)
  
  #Type <- as.vector(c(rep(0, nc), rep(1, nd)))
  #design <- model.matrix(~ Type)
  #fit <- lmFit(wholeDat, design)
  #fit <- eBayes(fit)
  #topTable(fit, coef = 1)
  #tT <- toptable(fit, genelist=rownames(data), number=nrow(data),adjust='fdr',p.value=1,lfc=log2(1))
  #resultFC <- tT[,c('ID','logFC','P.Value','adj.P.Val')]  
  
  ##########
  DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
  DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
  DEGenes <- foldChange[names(DEGenes)]
  
  peRes=pe(x=DEGenes,graphs=kpg,ref=names(foldChange),nboot=iterno)
  
  temp <- Summary(peRes)
  temp[,dim(temp)[2]+1] = rownames(temp)
  rownames(temp)<-kpn[rownames(temp)]
  colnames(temp)[dim(temp)[2]] = "pathway"
  temp
}

l <- list()
for (i in 1:loopNo) {
  a=sample(colnames(Samples),nc+nd,replace=FALSE)
  l[[i]] <- a
}

#library(multicore)
resNull <- mclapply(X=l,FUN=generateNull,mc.cores=20)

resNullORA=data.frame(row.names=names(kpn), Name=kpn,M=rep(NA,length(kpn)),std=rep(NA,length(kpn)))
for (i in 1:loopNo) {
  if (length(resNull[[i]]$pORA)==0) {resNullORA[,paste("pORA",i,sep="")]=NA}
  else {resNullORA[resNull[[i]]$pathway,paste("pORA",i,sep="")]=resNull[[i]]$pORA}
}
resNullORA$M=apply(resNullORA[4:(loopNo+3)],1,FUN=mean,na.rm=T)
resNullORA$std=apply(resNullORA[4:(loopNo+3)],1,FUN=sd,na.rm=T)

resNullIA=data.frame(row.names=names(kpn), Name=kpn,M=rep(NA,length(kpn)),std=rep(NA,length(kpn)))
for (i in 1:loopNo) {
  if (length(resNull[[i]]$pComb)==0) {resNullIA[,paste("pIA",i,sep="")]=NA}
  else {resNullIA[resNull[[i]]$pathway,paste("pIA",i,sep="")]=resNull[[i]]$pComb}
}
resNullIA$M=apply(resNullIA[4:(loopNo+3)],1,FUN=mean,na.rm=T)
resNullIA$std=apply(resNullIA[4:(loopNo+3)],1,FUN=sd,na.rm=T)

resNullPert=data.frame(row.names=names(kpn), Name=kpn,M=rep(NA,length(kpn)),std=rep(NA,length(kpn)))
for (i in 1:loopNo) {
  if (length(resNull[[i]]$pPert)==0) {resNullPert[,paste("pPert",i,sep="")]=NA}
  else {resNullPert[resNull[[i]]$pathway,paste("pPert",i,sep="")]=resNull[[i]]$pPert}
}
resNullPert$M=apply(resNullPert[4:(loopNo+3)],1,FUN=mean,na.rm=T)
resNullPert$std=apply(resNullPert[4:(loopNo+3)],1,FUN=sd,na.rm=T)

resNullAcc=data.frame(row.names=names(kpn), Name=kpn,M=rep(NA,length(kpn)),std=rep(NA,length(kpn)))
for (i in 1:loopNo) {
  if (length(resNull[[i]]$pAcc)==0) {resNullAcc[,paste("pAcc",i,sep="")]=NA}
  else {resNullAcc[resNull[[i]]$pathway,paste("pAcc",i,sep="")]=resNull[[i]]$pAcc}
}
resNullAcc$M=apply(resNullAcc[4:(loopNo+3)],1,FUN=mean,na.rm=T)
resNullAcc$std=apply(resNullAcc[4:(loopNo+3)],1,FUN=sd,na.rm=T)

file = paste("IA",nc,"vs",nd,"_400_",loopNo,".RData",sep="") 
save(resNullORA,resNullIA,resNullPert,resNullAcc,l,file=file)

PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullROntoTools.pdf"
pdf(file=PDFPath)  

#for (i in 1:length(resNullIA$Name))   
for (i in seq(150))  
{
  print(i)
  if(i == 41 || i == 42 ) {} 
  else {
    hist(as.numeric(resNullIA[i,-1]), breaks = 50, main = resNullIA[i,]$Name)  
  }
} 
dev.off() 



