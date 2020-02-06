path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
#dataSets=c("GSE1297","GSE28146","GSE5281_EC","GSE5281_HIP","GSE5281_MTG","GSE5281_PC", "GSE5281_SFG", "GSE5281_VCX")
dataSets=c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
           "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")

iterno=1000

library(PADOG)

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
# library("ROntoTools")
# kpg <- keggPathwayGraphs(organism)
# kpg <- setEdgeWeights(kpg)
# kpg <- setNodeWeights(kpg, defaultWeight = 1)
# kpn <- keggPathwayNames(organism)
load("KEGG65.150Pathways.RData")

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
	load(paste(path,"dataset/",dataset,"/",dataset,".RData",sep=""))
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

save(myGSARes,myGSAComb,dataSets,file="alldatasetsGSA.RData")

################ number of significant pathway ###########

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/alldatasetsGSA.RData")
signPW = vector(mode = "numeric", length = 0)
for (i in 1:length(myGSARes)) {
  signPW = c(signPW, length(which(myGSARes[[i]]$pGSA <0.05))/150)
}

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/GSA40Datasets.RData")
for (i in 2:(ncol(myGSARes)-3)) {
  signPW = c(signPW, length(which(myGSARes[,i]<0.05))/150)
}

mean(signPW)


