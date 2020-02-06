PVAL = 0.05
maxDE = 400

dataset = "GSE18842"
setwd("/Users/li11/myGit/pathwayReview/reviewDatasets/")
setwd("x:/myGit/pathwayReview/reviewDatasets/")


load(paste(dataset,"/",dataset,".RData",sep=""))   

data <- get(paste("gene_",dataset,sep=""))
group <- get(paste("group_",dataset,sep=""))

rownames(data) <- gsub("hsa:", "", rownames(data))
data <- as.matrix(data)

#Delete rows that data is constant for every sample
delRow <- sapply(seq(nrow(data)), function(x) { ifelse(sd(data[x,]) == 0, FALSE, TRUE) })
length(which(delRow))
data <- data[delRow,]

controlDat = data[,group$Group %in% 'c']
diseaseDat = data[,group$Group %in% 'd']

controlMean <- apply(controlDat,MARGIN=1,mean)
diseaseMean <- apply(diseaseDat,MARGIN=1,mean)

foldChange <- diseaseMean-controlMean

pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)

names(pvalues) <- names(foldChange)
DEGenes <- sort(abs(foldChange[pvalues<PVAL]),decreasing=TRUE)
DEGenes <- DEGenes[1:min(maxDE,length(DEGenes))]
DEGenes <- foldChange[names(DEGenes)]

##==============================
##  comparing dataset
##  non-small cell lung cancer
##==============================
#GSE19188
#GSE19188
#GSE19804
#GSE50627
#GSE6044

dataset = "GSE18842"
dataset = "GSE19188"
dataset = "GSE19804"
dataset = "GSE50627"
dataset = "GSE6044"

load(paste(dataset,"/",dataset,".RData",sep=""))   

data.GSE18842  <- get(paste("gene_",dataset,sep=""))
group.GSE18842 <- get(paste("group_",dataset,sep=""))
row.names(data.GSE18842)


data.GSE19188  <- get(paste("gene_",dataset,sep=""))
group.GSE19188 <- get(paste("group_",dataset,sep=""))

data.GSE19804  <- get(paste("gene_",dataset,sep=""))
group.GSE19804 <- get(paste("group_",dataset,sep=""))

data.GSE50627 <- get(paste("gene_",dataset,sep=""))
group.GSE50627<- get(paste("group_",dataset,sep=""))

data.GSE6044 <- get(paste("gene_",dataset,sep=""))
group.GSE6044<- get(paste("group_",dataset,sep=""))

length(which(row.names(data.GSE50627) %in% row.names(data.GSE6044)))

##  pancriatic cancer

dataset = "GSE28735"
load(paste(dataset,"/",dataset,".RData",sep=""))   

data.GSE28735  <- get(paste("gene_",dataset,sep=""))
group.GSE28735 <- get(paste("group_",dataset,sep=""))

dim(data.GSE28735)
length(which(row.names(data.GSE50627) %in% row.names(data.GSE28735)))



#GSE1418 Chronic Myeloid Leukemia 6 8 15618956 Bone marrow HG-Focus
dataset = "GSE1418"
data.GSE1418 <- get(paste("gene_",dataset,sep=""))
group.GSE1418 <- get(paste("group_",dataset,sep=""))
dim(data.GSE1418)
length(which(row.names(data.GSE6044) %in% row.names(data.GSE1418)))


annotation_GSE50627
annotation_GSE6044
annotation_GSE18842
