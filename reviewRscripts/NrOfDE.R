library(ggplot2)
path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/dataset/"
setwd(path)
NrDE = rep(NA, length(samePlat))
filename="datasetslist.txt"
dataSets = as.character(t(read.table(filename,header=F,sep="\t",stringsAsFactors=F)))
samePlat=as.list(dataSets)

for (i in 1:length(samePlat)) {
  dataset=samePlat[[i]]
  print(dataset)
  
  load(paste(path,dataset,"/",dataset,".RData",sep=""))      
  data <- get(paste("gene_",dataset,sep=""))
  group <- get(paste("group_",dataset,sep=""))
  
  controlDat = data[,group$Group %in% 'c']
  diseaseDat = data[,group$Group %in% 'd']
  
  controlMean <- apply(controlDat,MARGIN=1,mean)
  diseaseMean <- apply(diseaseDat,MARGIN=1,mean)
  
  foldChange <- diseaseMean-controlMean
  names(foldChange) <- paste("hsa:",names(foldChange),sep="")
  
  pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
  
  names(pvalues) <- names(foldChange)
  print("Classic")
  adj.pvalues <- p.adjust(pvalues, method = 'fdr')
  DEGenes <- foldChange[which(abs(foldChange[adj.pvalues<0.05]) > 1.5)]
  NrDE[i] <- length(DEGenes)
  
}

hist(NrDE, breaks=100)


##############################################################################################

load("/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/KEGG65.150Pathways.RData")
# load("/wsu/home/gd/gd03/gd0393/Pathway/KEGG65.150Pathways.RData")
kpg <- setEdgeWeights(kpg)
kpg <- setNodeWeights(kpg, defaultWeight = 1)

getnodes=function(x){
  return (x@nodes);}
path="/Volumes/MinhsBook/PathwayReview/dataset/"
setwd(path)

NrDE35 = rep(NA, 35)
samePlat=c("GSE781", "GSE1297","GSE3467","GSE3585","GSE3678","GSE4107","GSE5281_EC",
           "GSE5281_HIP","GSE5281_VCX","GSE6956_AA","GSE6956_C","GSE8671","GSE8762",
           "GSE9348","GSE9476","GSE14762","GSE15471","GSE16515","GSE18842","GSE19188",
           "GSE19728","GSE20153","GSE20291","GSE21354","GSE14924_CD4",
           "GSE14924_CD8","GSE16759","GSE19420","GSE20164","GSE23878","GSE24739_G0",
           "GSE24739_G1","GSE32676","GSE4183","GSE7305")
for (i in 1:length(samePlat)) {
  dataset = samePlat[i]
  print(dataset)
  load(paste0(dataset,"/",dataset, ".RData"))
  
  filename=paste(dataset,"/",dataset,"_annotation.txt",sep="")
  annotation=as.character((read.table(filename,header=F,sep="\t",stringsAsFactors=F))[1,1])
  keggnodes=lapply(kpg,getnodes)
  allGenes=unique(unlist(keggnodes))
  ano=paste(substring(annotation, 1, nchar(annotation)-3),"ENTREZID",sep="")
  require(annotation, character.only = TRUE)
  geneAno=as.data.frame(get(ano))
  geneAnno=geneAno[,2]
  names(geneAnno)=geneAno[,1]
  allProbes=unique(geneAno$probe_id)
  
  data <- get(paste("data_",dataset,sep=""))
  group <- get(paste("group_",dataset,sep=""))
  # data <- 2^data
  data <- data[rownames(data)%in%allProbes,]
  rownames(data) <- as.character(geneAnno[rownames(data)])
  data <- aggregate(data,by=list(rownames(data)),FUN=median)
  data <- data[paste("hsa:",data$Group.1,sep="")%in%allGenes,]
  rownames(data)=data$Group.1
  data <- data[,!colnames(data)%in%c("Group.1")]
  controlDat = data[,group$Group %in% 'c']
  diseaseDat = data[,group$Group %in% 'd']
  
  controlMean <- apply(controlDat,MARGIN=1,mean)
  diseaseMean <- apply(diseaseDat,MARGIN=1,mean)
  
  foldChange <- diseaseMean-controlMean
  names(foldChange) <- paste("hsa:",names(foldChange),sep="")
  
  pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
  
  names(pvalues) <- names(foldChange)
  adj.pvalues <- p.adjust(pvalues, method = 'fdr')
  #DEGenes <- foldChange[which(abs(foldChange[adj.pvalues<0.05]) > 1.5)]
  DEGenes <- intersect(names(which(adj.pvalues < 0.05)), names(which(abs(foldChange) > 1)))
  print(length(DEGenes))
  NrDE35[i] <- length(DEGenes)
}

NrDE75 = c(NrDE, NrDE35)
hist(NrDE75, breaks=100, main = "Histogram of number of DE genes", xlab = "Number of DE genes")
table(NrDE75)
save(NrDE, NrDE35, NrDE75, file = "NrDE75.RData")

load("NrDE75.RData")
ggplot(data = as.data.frame(NrDE), aes(NrDE)) + 
  geom_histogram(color="black", fill="darkgoldenrod1", bins = 40) +
  labs(x = "Number of DE genes", y = "Frequency") +
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold"))

kodatasets = c(2,0,6,0,1,0,20,0,27,1,0)
ggplot(data = as.data.frame(kodatasets), aes(kodatasets)) + 
  geom_histogram(color="black", fill="darkgoldenrod1") +
  labs(x = "Number of DE genes", y = "Frequency") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))






df <- data.frame(
  sex=factor(rep(c("F", "M"), each=200)),
  weight=round(c(rnorm(200, mean=55, sd=5), rnorm(200, mean=65, sd=5)))
)
head(df)

