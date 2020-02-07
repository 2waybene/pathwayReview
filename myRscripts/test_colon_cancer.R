##=====================================
##  test_colon_cancer.R
##====================================

#GSE4107	Colorectal cancer 	10	12	17317818	Colonic mucosa	HG-U133 Plus 2.0
#GSE8671	Colorectal cancer 	32	32	18171984	Colon      	HG-U133 Plus 2.0
#GSE9348	Colorectal cancer 	12	70	20143136	Colon      	HG-U133 Plus 2.0
#GSE23878	Colorectal cancer 	19	19	21281787	Colon      	HG-U133 Plus 2.0
#GSE4183	Colorectal cancer 	8	15	18776587	Colon      	HG-U133 Plus 2.0

#
BiocManager::install("hgu133a2.db")

library("hgu133plus2.db")
library("hgu95a.db")
library("hgu133a.db")
library("hgfocus.db")
library("hugene10sttranscriptcluster.db")
library(affyPLM)
library(affy)
library(simpleaffy)
library("GEOquery")
library("ROntoTools")
library("hgu133a2.db")

library(ggplot2)
library(ggthemes)
# library(psych)
library(data.table)



############### FIGURE 2: ranks and p-values of individual method ############### 
# library(vioplot)


setwd("x:/myGit/pathwayReview/")
path="x:/myGit/pathwayReview/"
filename="x:/myGit/pathwayReview/reviewDatasets/datasetslist.txt"

setwd("/Users/li11/myGit/pathwayReview/")
path="/Users/li11//myGit/pathwayReview/"
filename="/Users/li11//myGit/pathwayReview/reviewDatasets/datasetslist.txt"

rndataSets = as.character(t(read.table(filename,header=F,sep="\t",stringsAsFactors=F)))

methods <- c("FE", "WebGestalt", "GOstats", "KS", "WRS", "GSEA", "GSA", "PADOG", 
             "SPIA", "ROntoTools", "CePaGSA", "CePaORA", "PathNet") #change order
Type <- c("Non-TB", "Non-TB", "Non-TB", "Non-TB", "Non-TB", "Non-TB", "Non-TB", "Non-TB", 
          "TB", "TB", "TB", "TB", "TB")

############### FIGURE 2: ranks and p-values of individual method ############### 


dataSets=c("GSE4107", "GSE8671","GSE9348", "GSE23878","GSE4183")


#dataSet=c("GSE4107")
iterno=1000
PVAL <- 0.05
maxDE=400
load("reviewDatasets/KEGG65.150Pathways.RData")
kpg <- setEdgeWeights(kpg)
kpg <- setNodeWeights(kpg, defaultWeight = 1)

getnodes=function(x){
  return (x@nodes);}

keggnodes=lapply(kpg,FUN=function(x){return (x@nodes)})
allGenes=unique(unlist(keggnodes))

length(keggnodes)
#[1] 150
length(allGenes)
#4711

#mymapping=as.data.frame(org.Hs.egSYMBOL)
#mygeneSym=as.character(mymapping$symbol)
#names(mygeneSym)=mymapping$gene_id
#load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/mygeneSym.RData")
#load("/wsu/home/gd/gd03/gd0393/Pathway/mygeneSym.RData")

##==================
##  "ROntoTools",
##==================
myData <- list()
for(i in 1:(length(dataSets))) {
  
  dataset=dataSets[i]
  print(dataset)
  load(paste(path,"reviewDatasets/", dataset,"/",dataset,".RData",sep=""))
  
  group <- get(paste("group_",dataset,sep=""))
  annotation=get(paste("annotation_",dataset,sep=""))
  
  design=c("Not Paired")
  
  ano=paste(substring(annotation, 1, nchar(annotation)-3),"ENTREZID",sep="")
  require(annotation, character.only = TRUE)
  geneAno=as.data.frame(get(ano))
  geneAnno=geneAno[,2]
  names(geneAnno)=geneAno[,1]
  allProbes=unique(geneAno$probe_id)
  
  data <- 2^get(paste("data_",dataset,sep=""))
  data <- data[rownames(data)%in%allProbes,]
  rownames(data) <- as.character(geneAnno[rownames(data)])
  mydat <- aggregate(data,by=list(rownames(data)),FUN=median)
  mydat <- mydat[paste("hsa:",mydat$Group.1,sep="")%in%allGenes,]
  
  
  controlDat = mydat[,rownames(group[group$Group=="c",])]
  rownames(controlDat) = as.character(mydat$Group.1)
  
  diseaseDat = mydat[,rownames(group[group$Group=="d",])]
  rownames(diseaseDat) = as.character(mydat$Group.1)
  
  
  controlMean <- log(apply(controlDat,MARGIN=1,mean),base=2)
  diseaseMean <- log(apply(diseaseDat,MARGIN=1,mean),base=2)
  
  foldChange <- diseaseMean-controlMean
  names(foldChange) <- paste("hsa:",names(foldChange),sep="")
  
  getTtestPvalue <- function(x, controlDat, diseaseDat) {  
    if(var(t(controlDat[x,])) ==0  | var (t(diseaseDat[x,])) == 0 ) {return(1)}
    else{
      return ( t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
    }
  }
  
  #if paired, needs some attention
  pvalues <- sapply(seq(nrow(controlDat)), function (x) {getTtestPvalue(x, controlDat, diseaseDat)}) 
  names(pvalues) <- names(foldChange)
  
  #DEGenes <- foldChange[(abs(foldChange)>FC) & (pvalues<PVAL)]
  DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
  #maxDE=min(floor(length(allGenes)/10),length(DEGenes))
  DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
  DEGenes <- foldChange[names(DEGenes)]
  
  
  length(DEGenes)
  
  peRes=pe(x=DEGenes,graphs=kpg,ref=names(foldChange),nboot=iterno)
  temp <-Summary(peRes, pathNames = kpn, totalAcc = FALSE, totalPert = FALSE,
               pAcc = FALSE, pORA = FALSE, comb.pv = NULL, order.by = "pPert")
  
#  temp <- Summary(peRes)
 # temp <- cbind.data.frame(temp,kpn[rownames(temp)])
  #colnames(temp) <- c(colnames(Summary(peRes)), "pathway")
  
  # XX <- rep(NA,length(kpn))
  # temp1=data.frame(row.names=names(kpn), Name=kpn, pPert=XX, 
  #                  pPert.adj=XX, )
  
  #temp <- cbind.data.frame(temp$pathway, temp$pComb)
  # rownames(temp) <- rownames(Summary(peRes))
  # temp = na.omit(temp)
  # head(temp)
  myData[[i]] <- temp
}

names(myData) <- dataSets
save(myData,file=paste(path, "researchDatasets/colonCancer_ROntoTools.RData",sep=""))

## 
##==================
##  "FE",
##==================
myData <- list()
for(i in 1:(length(dataSets))) {
  
  dataset=dataSets[i]
  print(dataset)
  load(paste(path,"reviewDatasets/", dataset,"/",dataset,".RData",sep=""))

  group <- get(paste("group_",dataset,sep=""))
  annotation=get(paste("annotation_",dataset,sep=""))

  design=c("Not Paired")
  
  ano=paste(substring(annotation, 1, nchar(annotation)-3),"ENTREZID",sep="")
  require(annotation, character.only = TRUE)
  geneAno=as.data.frame(get(ano))
  geneAnno=geneAno[,2]
  names(geneAnno)=geneAno[,1]
  allProbes=unique(geneAno$probe_id)
  
  data <- 2^get(paste("data_",dataset,sep=""))
  data <- data[rownames(data)%in%allProbes,]
  rownames(data) <- as.character(geneAnno[rownames(data)])
  mydat <- aggregate(data,by=list(rownames(data)),FUN=median)
  mydat <- mydat[paste("hsa:",mydat$Group.1,sep="")%in%allGenes,]
  
  
  controlDat = mydat[,rownames(group[group$Group=="c",])]
  rownames(controlDat) = as.character(mydat$Group.1)
  
  diseaseDat = mydat[,rownames(group[group$Group=="d",])]
  rownames(diseaseDat) = as.character(mydat$Group.1)
  
  
  controlMean <- log(apply(controlDat,MARGIN=1,mean),base=2)
  diseaseMean <- log(apply(diseaseDat,MARGIN=1,mean),base=2)
  
  foldChange <- diseaseMean-controlMean
  names(foldChange) <- paste("hsa:",names(foldChange),sep="")
  
  getTtestPvalue <- function(x, controlDat, diseaseDat) {  
    if(AllEqual(t(controlDat[x,])) | AllEqual(t(diseaseDat[x,]))) {return(1)}
    else{
      return ( t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
    }
  }

  #if paired, needs some attention
  pvalues <- sapply(seq(nrow(controlDat)), function (x) {getTtestPvalue(x, controlDat, diseaseDat)}) 
  names(pvalues) <- names(foldChange)
  
  #DEGenes <- foldChange[(abs(foldChange)>FC) & (pvalues<PVAL)]
  DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
  #maxDE=min(floor(length(allGenes)/10),length(DEGenes))
  DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
  DEGenes <- foldChange[names(DEGenes)]
  
  temp=data.frame(row.names=names(kpn), Name=kpn,pvalue=rep(NA,length(kpn)),rank=rep(NA,length(kpn)))
  for (pw in 1:length(kpn)) {
    #print(kpn[pw])
    DEinPW = length(intersect(keggnodes[[pw]],names(DEGenes)))
    #print(DEinPW)
    DEnotInPW = length(DEGenes) - DEinPW
    notDEinPW = length(keggnodes[[pw]]) - DEinPW
    notDEnotinPW = length(pvalues) - length(DEGenes) - length(keggnodes[[pw]]) + DEinPW
    table <- matrix(c(DEinPW, DEnotInPW, notDEinPW, notDEnotinPW), nrow = 2, dimnames = list(Genes = c("in PW", "not in PW"),Pathway = c("DE", "not DE")))
    fisher <- fisher.test(table, alternative = "greater")
    #fisher$p.value
    temp[pw,"pvalue"] <- fisher$p.value
  }
  temp$rank <- rank(temp[,2],ties.method="min")
  myData[[i]] <- temp
}

names(myData) <- dataSets
save(myData,file=paste(path, "researchDatasets/colonCancer_FE.RData",sep=""))

##=================
##  GSA
##=================
library("PADOG")
iterno=1000



fishersMethod = function(x) {
  pchisq(-2 * sum(log(x)), df=2*length(x), lower=FALSE)
}

IrwinHall = function(x,n) {
  1/factorial(n) * sum(sapply(0:n, function(k) (-1)^k * choose(n,k) * ifelse(x > k,x-k,0)^(n)))
}

additiveMethod = function(x) {
  n = length(x)
  if (n <= 30) {
    IrwinHall(sum(x),n)
  } else {
    pnorm(sum(x),n/2,sqrt(n/12),lower=TRUE)
  }
}

pValCombMethods=list(Fisher=fishersMethod,Edgington=additiveMethod)


##=======================
##  GSA function
##=======================


gsaF = function(dat.m, ano, annotation, design="Not Paired", mygslist, gs.names, minsize=3) {
  require(GSA)
  
  if (design=="") {
    design="Not Paired"
  }
  group = ano$Group
  block = ano$Block
  esetm = dat.m
  if (!is.null(annotation)) {
    aT1 = filteranot(esetm = esetm, group = group, paired = (design == "Paired"), block, annotation)
    aT1 <- aT1[aT1$ENTREZID %in% (unlist(mygslist)),]
    esetm = esetm[rownames(esetm) %in% aT1$ID, ]
    rownames(esetm) <- aT1$ENTREZID[match(rownames(esetm), aT1$ID)]
  }
  nc = table(ano$Group)["c"]
  nd = table(ano$Group)["d"]
  if (design == "Not Paired") {
    yy = c(rep(1, nc), rep(2, nd))
  } else {
    block = as.numeric(factor(ano$Block))
    block[duplicated(block)] <- (-block[duplicated(block)])
    yy = block
  }
  resgsa = GSA(x = esetm, y = yy, genesets = mygslist, 
               genenames = rownames(esetm), method = "maxmean", 
               resp.type = ifelse(design == "Not Paired", "Two class unpaired", 
                                  "Two class paired"), censoring.status = NULL, 
               random.seed = 1, knn.neighbors = 10, s0 = NULL, s0.perc = NULL, 
               minsize = minsize, maxsize = 1000, restand = TRUE, 
               restand.basis = c("catalog", "data"), nperms = iterno)
  res = data.frame(Name=gs.names[names(mygslist)],ID = names(mygslist), pGSA = 2 * apply(cbind(resgsa$pvalues.lo, resgsa$pvalues.hi), 1, min), stringsAsFactors = FALSE)
  res <- na.omit(res)
  res = res[order(res$pGSA), ]
  rownames(res) <- res$ID
  res
}


organism="hsa"
gslist=lapply(kpg,FUN=function(x){return (x@nodes);})
gs.names=kpn[names(gslist)]
names(gs.names)=substr(names(gs.names), nchar("path:hsa")+1, 1000000L)
names(gslist)=substr(names(gslist), nchar("path:hsa")+1, 1000000L)

for (i in 1:length(gslist)) {
  gslist[[i]] = substr(gslist[[i]], nchar("hsa:")+1, 1000000L)
}


myGSARes=list()
for(i in 1:(length(dataSets))) {
  dataset=dataSets[i]
  load(paste(path,"reviewDatasets/", dataset,"/",dataset,".RData",sep=""))
  
  annotation=get(paste("annotation_",dataset,sep=""))
  design=get(paste("design_",dataset,sep=""))
  data <- get(paste("data_",dataset,sep=""))
  group <- get(paste("group_",dataset,sep=""))
  group <- group[order(group$Group),]
  data <- data[,rownames(group)]
  
  res <- gsaF (dat.m=data,
               ano=group, 
               annotation=annotation, 
               design=design, 
               mygslist=gslist, 
               gs.names=gs.names,
               minsize=3)
  myGSARes[[i]]=res
}


r=rownames(myGSARes[[1]]);
for (i in 2:length(myGSARes)) {
  r=r[r%in%rownames(myGSARes[[i]])]
}

for (i in 1:length(myGSARes)) {
  myGSARes[[i]]=myGSARes[[i]][rownames(myGSARes[[i]])%in%r,]
}

for (i in 1:length(myGSARes)){
  myGSARes[[i]][,"rGSA"]=seq(1,length(rownames( myGSARes[[i]])),1)
}

myGSAComb=data.frame(matrix(NA,nrow=length(rownames(myGSARes[[1]])),ncol=1+2*length(dataSets)+3*length(pValCombMethods)))
rownames(myGSAComb)=rownames(myGSARes[[1]])
colnames(myGSAComb)=c("Name",paste("pGSA",dataSets, sep="_"),paste("rGSA",dataSets, sep="_"),paste("p",names(pValCombMethods),sep=""), paste("p",names(pValCombMethods),".fdr",sep=""),paste("r",names(pValCombMethods),sep=""))
myGSAComb[,"Name"] = myGSARes[[1]][,"Name"]
for (i in 1:length(dataSets)) {
  myGSAComb[,paste("pGSA",dataSets[i],sep="_")] = myGSARes[[i]][rownames(myGSAComb),"pGSA"]
  myGSAComb[,paste("rGSA",dataSets[i],sep="_")] = myGSARes[[i]][rownames(myGSAComb),"rGSA"]
}

for (meth in names(pValCombMethods)){
  myGSAComb[,paste("p",meth,sep="")] = apply(myGSAComb[,2:(1+length(myGSARes))],1,pValCombMethods[[meth]])
  myGSAComb[,paste("p",meth,".fdr",sep="")] = p.adjust(myGSAComb[,paste("p",meth,sep="")],method="fdr")
  myGSAComb[order(myGSAComb[,paste("p",meth,sep="")]),paste("r",meth,sep="")] = seq(1,length(rownames(myGSAComb)),1)
}

myGSAComb=myGSAComb[order(myGSAComb$pEdgington),]

save(myGSARes,myGSAComb,dataSets,file=paste(path, "researchDatasets/colonCancer_GSA.RData",sep=""))

################ number of significant pathway ###########
load("x:/myGit/pathwayReview/researchDatasets/colonCancer_GSA.RData")
signPW = vector(mode = "numeric", length = 0)
for (i in 1:length(myGSARes)) {
  signPW = c(signPW, length(which(myGSARes[[i]]$pGSA <0.05))/150)
}

mean(signPW)

##=======================
##  GSEA function
##=======================
library(GSEA)

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

classFileCreation <- function (group, datasetdir="gsea_stored_datasets")
{
  tempCLS <- paste0 ( datasetdir, "/tempCLS.cls")
  sink (tempCLS)
  total <- length(group$Group)
  groups <- length(levels(as.factor(group$Group)))
  cat (total, " ", groups, " 1")
  cat ("\n")
  cat ("# ", levels(as.factor(group$Group)))
  cat ("\n")
  cat( t(group$Group))
  sink()
  return(tempCLS)
}

library(ROntoTools)
load("reviewDatasets/KEGG65.150Pathways.RData")
keggnodes=lapply(kpg,FUN=function(x){return (x@nodes)})
allGenes=unique(unlist(keggnodes));allGenes=substr(allGenes, nchar("hsa:")+1, 1000000L)

myGSEAComb=data.frame(matrix(NA, nrow=length(kpn), ncol=1+2*length(dataSets)+3*length(pValCombMethods)))
colnames(myGSEAComb)=c("Name",paste("pGSEA",dataSets, sep="_"),paste("rGSEA",dataSets, sep="_"),paste("p",names(pValCombMethods),sep=""), paste("p",names(pValCombMethods),".fdr",sep=""),paste("r",names(pValCombMethods),sep=""))
rownames(myGSEAComb)=names(kpn)
myGSEAComb[,"Name"]=kpn


myData <- list()
for (i in 1:(length(dataSets))) {
  dataset=dataSets[i]
  load(paste(path,"reviewDatasets/", dataset,"/",dataset,".RData",sep=""))
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

save(myData,file=paste(path, "researchDatasets/colonCancer_GSEA.RData",sep=""))
