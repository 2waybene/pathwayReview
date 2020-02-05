#setwd("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/");
#path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/"
setwd("/wsu/home/gd/gd03/gd0393/Pathway/dataset/");
path="/wsu/home/gd/gd03/gd0393/Pathway/dataset/"
dataSets=c("GSE3467","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE8671","GSE8762",
           "GSE9348","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")

iterno=1000

library(PADOG)

gsaF = function(dat.m, ano, annotation, design="Not Paired", mygslist, minsize=3) {
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
library("ROntoTools")
kpg <- keggPathwayGraphs(organism)
kpg <- setEdgeWeights(kpg)
kpg <- setNodeWeights(kpg, defaultWeight = 1)
kpn <- keggPathwayNames(organism)

gslist=lapply(kpg,FUN=function(x){return (x@nodes);})
gs.names=kpn[names(gslist)]
names(gs.names)=substr(names(gs.names), nchar("path:hsa")+1, 1000000L)
names(gslist)=substr(names(gslist), nchar("path:hsa")+1, 1000000L)

for (i in 1:length(gslist)) {
  gslist[[i]] = substr(gslist[[i]], nchar("hsa:")+1, 1000000L)
}


setwd("/Users/li11/myGit/pathwayReview")
path ="reviewDatasets/"
#generate the null distribution
samePlat=c("GSE3467","GSE3678","GSE4107","GSE5281EC",
           "GSE5281HIP","GSE5281VCX","GSE8671","GSE8762",
           "GSE9348","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE21354","GSE14924CD4",
           "GSE14924CD8","GSE16759","GSE19420","GSE23878","GSE24739G0",
           "GSE24739G1","GSE32676","GSE4183","GSE7305")
annotation="hgu133plus2.db"
ctrSamples <- NULL
for (i in 1:length(samePlat)) {
  dataset=samePlat[i]
  load(paste(path,dataset,"/",dataset,".RData",sep=""))      
  data <- get(paste("data_",dataset,sep=""))
  group <- get(paste("group_",dataset,sep=""))
  cn = rownames(group[group$Group=="c",])
  ctrSamples <- cbind(ctrSamples,data[,cn])
}


#args <- commandArgs(trailingOnly = TRUE)
#args <- as.numeric(args)
#nc=args[1]
#nd=args[2]
#loopNo=args[3]

nc = 15
nd = 15
loopNo = 5


generateNull = function (a) {
  group <- data.frame(row.names=a,Sample=a,Group=as.factor(c(rep("c",nc),rep("d",nd))))
  data=ctrSamples[,a]
  res <- gsaF (dat.m=ctrSamples[,a],
               ano=data.frame(row.names=a,Sample=a,Group=as.factor(c(rep("c",nc),rep("d",nd)))), 
               annotation=annotation, 
               design="Not Paired", 
               mygslist=gslist, 
               minsize=3)
}

l <- list()
for (i in 1:loopNo) {
  a=sample(colnames(ctrSamples),nc+nd,replace=FALSE)
  l[[i]] <- a
}

library(multicore)
resNull <- mclapply(X=l,FUN=generateNull,mc.cores=10)
resNull <- lapply(X=l,FUN=generateNull)
r=rownames(resNull[[1]]);
for (i in 2:length(resNull)) {
  r=r[r%in%rownames(resNull[[i]])]
}
for (i in 1:length(resNull)) {resNull[[i]]=resNull[[i]][r,]} 
resNullAll=data.frame(row.names=r, Name=resNull[[1]][r,"Name"],M=rep(NA,length(r)),std=rep(NA,length(r)))
for (i in 1:loopNo) {resNullAll[,paste("pGSA",i,sep="")]=resNull[[i]][rownames(resNullAll),"pGSA"]}
resNullAll[,"M"]=apply(resNullAll[4:(loopNo+3)],1,FUN=mean)
resNullAll[,"std"]=apply(resNullAll[4:(loopNo+3)],1,FUN=sd)
file = paste("GSA",nc,"vs",nd,"_",loopNo,".RData",sep="") 
save(resNullAll,l,file=file)
