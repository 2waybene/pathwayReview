rm(list=ls())
setwd("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/")
path= "/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
source("gsea.files.creation.R")
source("GSEA.1.0.R", verbose=FALSE, max.deparse.length=9999)

source("DANUBE_functions.R")

getGeneExpression = function(data, annotation, allGenes=NULL) {
  ano=paste(substring(annotation, 1, nchar(annotation)-3),"ENTREZID",sep="")
  require(annotation, character.only = TRUE)
  geneAno=as.data.frame(get(ano))
  if (!is.null(allGenes)) {geneAno=geneAno[geneAno$gene_id%in%allGenes,]}
  geneAnno=geneAno[,2]
  names(geneAnno)=geneAno[,1]
  allProbes=unique(geneAno$probe_id)
  
  data = 2^data
  data <- data[rownames(data)%in%allProbes,]
  rownames(data) <- as.character(geneAnno[rownames(data)])
  
  mydat <- aggregate(data,by=list(rownames(data)),FUN=median)
  rownames(mydat)=mydat$Group.1
  
  mydat <- log(mydat[,-1],2)
}

library(ROntoTools)
load("KEGG65.150Pathways.RData")
keggnodes=lapply(kpg,FUN=function(x){return (x@nodes)})
allGenes=unique(unlist(keggnodes));allGenes=substr(allGenes, nchar("hsa:")+1, 1000000L)

# dataSets=c("GSE28146","GSE5281_EC","GSE5281_HIP","GSE5281_MTG","GSE5281_PC", "GSE5281_SFG", "GSE5281_VCX")
# dataSets =c("GSE781", "GSE1297")
dataSets=c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
           "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")

myGSEAComb=data.frame(matrix(NA, nrow=length(kpn), ncol=1+2*length(dataSets)+3*length(pValCombMethods)))
colnames(myGSEAComb)=c("Name",paste("pGSEA",dataSets, sep="_"),paste("rGSEA",dataSets, sep="_"),paste("p",names(pValCombMethods),sep=""), paste("p",names(pValCombMethods),".fdr",sep=""),paste("r",names(pValCombMethods),sep=""))
rownames(myGSEAComb)=names(kpn)
myGSEAComb[,"Name"]=kpn

myData <- list()
for (i in 1:(length(dataSets))) {
  dataset=dataSets[i]
  load(paste(path, "dataset/",dataset,"/",dataset, ".RData", sep=""))
  group <- get(paste("group_",dataset,sep=""))
  annotation=get(paste("annotation_",dataset,sep=""))
  data=get(paste("data_",dataset,sep=""))
  
  data=getGeneExpression(data, annotation, allGenes)
  rownames(data)=paste("hsa:", rownames(data), sep="") 
  data=data[,group$Sample]
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
  res=res[order(res$"NOM p-val"),]
  res$pval.fdr=p.adjust(as.numeric(paste(res$"NOM p-val")), method="fdr")
  res$rank=seq(nrow(res))
  
  # myGSEAComb[rownames(res),paste("pGSEA_",dataset, sep="")]=as.numeric(paste(res$"NOM p-val"))
  # myGSEAComb[rownames(res),paste("rGSEA_",dataset, sep="")]=res$rank
  
  temp <- res[,c(3,6,15)]
  colnames(temp) = c("Pathway", "pvalue", "rank")
  myData[[i]] <- temp
  # myGSEARes[[length(myGSEARes)+1]]=res;
}
names(myData) = dataSets

for (meth in names(pValCombMethods)){
  myGSEAComb[,paste("p",meth,sep="")] = apply(myGSEAComb[,2:(1+length(dataSets))],1,pValCombMethods[[meth]])
  myGSEAComb[,paste("p",meth,".fdr",sep="")] = p.adjust(myGSEAComb[,paste("p",meth,sep="")],method="fdr")
  myGSEAComb[order(myGSEAComb[,paste("p",meth,sep="")]),paste("r",meth,sep="")] = seq(1,length(rownames(myGSEAComb)),1)
}
myGSEAComb=myGSEAComb[order(rownames(myGSEAComb)),]
myGSEAComb=myGSEAComb[order(myGSEAComb[,"pEdgington"]),]

file="alldatasetsGSEA.RData"
save(myData,file=file)

################ number of significant pathway ###########

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsGSEA.RData")
signPW = vector(mode = "numeric", length = 0)
for (i in 1:length(myData)) {
  temp <- as.vector(myData[[i]]$pvalue)
  signPW = c(signPW, length(which(temp <0.05))/150)
}

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/GSEA40Datasets.RData")
for (i in 1:length(myData)) {
  temp <- as.vector(myData[[i]]$pvalue)
  signPW = c(signPW, length(which(temp<0.05))/150)
}

mean(signPW)

