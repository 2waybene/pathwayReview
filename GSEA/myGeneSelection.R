PVAL = 0.05
maxDE = 400

dataset = "GSE18842"
setwd("x:/myGit/pathwayReview/reviewDatasets/")
load(paste(dataset,"/",dataset,".RData",sep="")) 
data <- get(paste("gene_",dataset,sep=""))
group <- get(paste("group_",dataset,sep=""))

  



rownames(data) <- gsub("hsa:", "", rownames(data))
data <- as.matrix(data)

#Delete rows that data is constant for every sample
delRow <- sapply(seq(nrow(data)), function(x) { ifelse(sd(data[x,]) == 0, TRUE, FALSE) })
if (length(which(delRow)) > 0)
{
  data <- data[-c(which(delRow)),]
}

#group <- get(paste("group_",dataset,sep=""))

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


