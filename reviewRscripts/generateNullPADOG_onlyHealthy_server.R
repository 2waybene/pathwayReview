#setwd("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/");
#path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/"
setwd("/wsu/home/gd/gd03/gd0393/Pathway/dataset/");
path="/wsu/home/gd/gd03/gd0393/Pathway/dataset/"
#dataSets=c("GSE1297","GSE28146","GSE5281_EC","GSE5281_HIP","GSE5281_MTG","GSE5281_PC", "GSE5281_SFG", "GSE5281_VCX")
dataSets=c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
           "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")
iterno=1000
targetGeneSets="05010"

library(PADOG)


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


#generate the null distribution
#samePlat=c("GSE28146","GSE5281_EC","GSE5281_HIP","GSE5281_MTG","GSE5281_PC", "GSE5281_SFG", "GSE5281_VCX")
samePlat = c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
             "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
             "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
             "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
             "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
             "GSE24739_G1","GSE32676","GSE4183","GSE7305")
annotation="hgu133plus2.db"
ctrSamples <- NULL
intersectGenes <- NULL
for (i in 1:length(samePlat)) {
  dataset=samePlat[i]
  print(dataset)
  load(paste(path,dataset,"/",dataset,".RData",sep=""))      
  data <- get(paste("data_",dataset,sep=""))
  group <- get(paste("group_",dataset,sep=""))
  cn = rownames(group[group$Group=="c",])
  
  if (is.null(intersectGenes)) {
    intersectGenes = rownames(data)
  } else {
    intersectGenes = intersect(intersectGenes, rownames(data))
  }
  if (is.null(ctrSamples)) {
    ctrSamples = data
  } else {
    ctrSamples <- cbind(ctrSamples[intersectGenes,],data[intersectGenes,])
  }
}


args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
#nc=args[1]
#nd=args[2]
#loopNo=args[3]
nc = 15
nd = 15
loopNo = 5000

generateNull = function (a) {
  res = try(padog(
		esetm=ctrSamples[,a],
		group=as.factor(c(rep("c",nc),rep("d",nd))),
		paired=FALSE,
		block=NULL,
		targetgs=targetGeneSets,
		annotation=annotation,
		gslist=gslist,
		gs.names=gs.names,
		organism="hsa",
		verbose=TRUE,
		Nmin=3,
		NI=iterno,
		plots=FALSE), silent=TRUE)
  res
}

l <- list()
for (i in 1:loopNo) {
  a=sample(colnames(ctrSamples),nc+nd,replace=FALSE)
  l[[i]] <- a
}

#library(multicore)
resNull <- mclapply(X=l,FUN=generateNull, mc.cores=20)
#resNull <- lapply(X=l,FUN=generateNull)

r=rownames(resNull[[1]]);
for (i in 2:length(resNull)) {
  r=union(r,rownames(resNull[[i]]))
}

resNullAll=data.frame(row.names=r, Name=gs.names[r], M=rep(NA,length(r)),std=rep(NA,length(r)))
for (i in 1:loopNo) {
	if (class(resNull[[i]])=="try-error") {resNullAll[,paste("pPADOG",i,sep="")]=NA}
	else {resNullAll[,paste("pPADOG",i,sep="")]=resNull[[i]][rownames(resNullAll),"Ppadog"]}
}

resNullAll$M=apply(resNullAll[4:(loopNo+3)],1,FUN=mean,na.rm=T)
resNullAll$std=apply(resNullAll[4:(loopNo+3)],1,FUN=sd,na.rm=T)
file = paste("PADOG",nc,"vs",nd,"_",loopNo,".RData",sep="") 
save(resNullAll,l,file=file)

