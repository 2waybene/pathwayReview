PVAL = 0.05
maxDE = 400

dataset = "GSE32676"
setwd("/Users/minhnguyen/75Datasets_PathwayReview/")
load(paste(dataset,"/",dataset,".RData",sep=""))      
data <- get(paste("gene_",dataset,sep=""))
group <- get(paste("group_",dataset,sep=""))

rownames(data) <- gsub("hsa:", "", rownames(data))
data <- as.matrix(data)

#Delete rows that data is constant for every sample
delRow <- sapply(seq(nrow(data)), function(x) { ifelse(sd(data[x,]) == 0, FALSE, TRUE) })
data <- data[delRow,]

group <- get(paste("group_",dataset,sep=""))

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


