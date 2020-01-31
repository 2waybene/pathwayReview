#source("http://bioconductor.org/biocLite.R")
#biocLite("CePa")
setwd("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/")
path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
# setwd("/wsu/home/gd/gd03/gd0393/Pathway/dataset/");
# path="/wsu/home/gd/gd03/gd0393/Pathway/dataset/"
library("ROntoTools")
source("gsea.files.creation.R")
source("GSEA.1.0.R", verbose=FALSE, max.deparse.length=9999)

source("DANUBE_functions.R")

iterno=2000
PVAL <- 0.05
maxDE=400

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
# samePlat=c("GSE7305")

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/mygeneSym.RData")
# load("/wsu/home/gd/gd03/gd0393/Pathway/mygeneSym.RData")

Samples <- NULL
intersectGenes <- NULL
for (i in 1:length(samePlat)) {
  dataset=samePlat[i]
  print(dataset)
  
  filename=paste(path,"dataset/",dataset,"/",dataset,"_annotation.txt",sep="")
  # filename=paste(path,dataset,"/",dataset,"_annotation.txt",sep="")
  
  annotation=as.character((read.table(filename,header=F,sep="\t",stringsAsFactors=F))[1,1])
  keggnodes=lapply(kpg,getnodes)
  allGenes=unique(unlist(keggnodes))
  ano=paste(substring(annotation, 1, nchar(annotation)-3),"ENTREZID",sep="")
  require(annotation, character.only = TRUE)
  geneAno=as.data.frame(get(ano))
  geneAnno=geneAno[,2]
  names(geneAnno)=geneAno[,1]
  allProbes=unique(geneAno$probe_id)
  load(paste(path,"dataset/",dataset,"/",dataset,".RData",sep=""))  
  # load(paste(path,dataset,"/",dataset,".RData",sep=""))      
  data <- get(paste("data_",dataset,sep=""))
  
  group <- get(paste("group_",dataset,sep=""))
  cn = rownames(group[group$Group=="c",])
  data <- data[,cn]
  
  # data <- 2^data
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
# nc=5
# nd=5
# loopNo=2
nc=15
nd=15
#loopNo=5000
loopNo=2000

generateNull = function (a) {
  data= Samples[,a]
  rownames(data)=paste("hsa:", rownames(data), sep="") 
  # data=data[,group$Sample]
  group <- data.frame(Sample = a, Group = c(rep("c", nc), rep("d", nd)))
  rownames(group) = a
  clsfile <- classFileCreation(group, datasetdir="gsea_stored_datasets")
  
  gseares <- GSEA(input.ds = data,
                  input.cls = clsfile,
                  gs.db = "gene_sets_db/tinkegg_ROT.gmt",
                  #output.directory = "gsea_stored_results/",
                  output.directory = "",
                  doc.string= "kegg",
                  non.interactive.run   = F,               
                  reshuffling.type      = "sample.labels", 
                  nperm                 = 1000,
                  weighted.score.type   =  1,            
                  nom.p.val.threshold   = -1,            
                  fwer.p.val.threshold  = -1,            
                  fdr.q.val.threshold   = 0.25,          
                  topgs                 = 20,            
                  adjust.FDR.q.val      = F,             
                  gs.size.threshold.min = 1,             
                  gs.size.threshold.max = 500,           
                  reverse.sign          = F,             
                  preproc.type          = 0,             
                  random.seed           = 3338,          
                  perm.type             = 0,             
                  fraction              = 1.0,           
                  replace               = F,             
                  save.intermediate.results = F,         
                  OLD.GSEA              = F,             
                  use.fast.enrichment.routine = T)
  
  
  rownames(gseares$report1)=gseares$report1$GS
  rownames(gseares$report2)=gseares$report2$GS
  
  res=rbind(gseares$report1, gseares$report2)
  res=res[order(res$GS),]
  res = res[,c("SOURCE", "NOM p-val")]
  colnames(res) = c("Name", "pvalue")
  res
}

l <- list()
for (i in 1:loopNo) {
  a=sample(colnames(Samples),nc+nd,replace=F)
  l[[i]] <- a
}

#library(multicore)
# resNull <- mclapply(X=l,FUN=generateNull,mc.cores=10)
resNull <- lapply(X=l,FUN=generateNull)



#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/CePaORA15vs15_400_2000 (original).RData")

noLoop <- length(resNull)
pathName <- resNull[[1]][,"Name"]

resNullGSEA=data.frame(row.names=pathName, Name=pathName,M=rep(NA,length(pathName)),std=rep(NA,length(pathName)))
for (i in 1:noLoop) {
  if (length(resNull[[i]][,"pvalue"])==0) {resNullGSEA[,paste("p_",i,sep="")]=NA}
  else {resNullGSEA[as.vector(resNull[[i]][,"Name"]),paste("p_",i,sep="")]=resNull[[i]][,"pvalue"]}
}
resNullGSEA$M=apply(resNullGSEA[4:(noLoop+3)],1,FUN=mean,na.rm=T)
resNullGSEA$std=apply(resNullGSEA[4:(noLoop+3)],1,FUN=sd,na.rm=T)

file = paste0("GSEA",nc,"vs",nd,"_400_",loopNo,".RData",sep="") 
save(resNullGSEA,l,file=file)


# load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/GSEA15vs15_400_2000.RData")
# PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullGSEA.pdf"
# pathName <- resNullGSEA[,"Name"]
# pdf(file=PDFPath)
# 
# for (i in seq(length(pathName)))  {
#   print(i)
#   hist(as.numeric(resNullGSEA[i,-c(1,2,3)]), breaks = 50, main = resNullGSEA[i,]$Name)
# }
# dev.off()
# 
# PDFPath = "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/generateNullGSEA_color.pdf"
# 
# pdf(file=PDFPath)
# 
# for (i in seq(length(pathName)))  {
#   print(i)
#   hist(as.numeric(resNullGSEA[i,-c(1,2,3)]), breaks = 50, main = resNullGSEA[i,]$Name, col = "blue")
# }
# dev.off()
# 
# 
