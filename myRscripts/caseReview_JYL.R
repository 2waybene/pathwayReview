# library(vioplot)
library(ggplot2)
library(ggthemes)
# library(psych)
library(data.table)

############### FIGURE 2: ranks and p-values of individual method ############### 
# library(vioplot)
setwd("x:/myGit/pathwayReview/")
path="x:/myGit/pathwayReview/"

dataFileInfo <- read.csv ("helpDoc/geo_75_downloads.csv", header = TRUE)
levels(dataFileInfo$Platform)

dataFileInfo[ order(dataFileInfo$Platform),]

case <- c()

case[1] <- "GSE6044"  ## HG-Focus
case[2] <- "GSE14762" ## HG-U133 Plus 2.0
case[3] <- "GSE781"   ## HG-U133A
case[4] <- "GSE39825" ## HG U95Av2
case[5] <- "GSE79962" ## HuGene-10st


data <- c()
for (i in 1:length(case))
{
  load(paste  ("reviewDatasets","/", case[i],"/",case[i] ,".RData",sep=""))
  data[i] <- list (dataFile = get(paste("gene_",case[i] ,sep="")))
}


str(data)

dim(data_GSE6044)
#[1] 8793   36
dim(gene_GSE6044)
#[1] 3101   36

dim(data_GSE14762)
#[1] 42323    22
dim(gene_GSE14762)
#[1] 4154   22

dim(data_GSE781)
#[1] 19855    17
dim(gene_GSE781)
#[1] 3671   17

dim(data_GSE39825)
#[1] 12625    10
dim(gene_GSE39825)
#[1] 3061   10

dim(data_GSE79962)
#[1] 32321    20
dim(gene_GSE79962)
#[1] 4372   20