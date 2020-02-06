
rm(list=ls())
library(GSA)



##======================
##  a single sample
##======================


dataset = "GSE18842"
setwd("x:/myGit/pathwayReview/reviewDatasets/")
load(paste(dataset,"/",dataset,".RData",sep=""))   

ls()
#[1] "annotation_GSE18842" "data_GSE18842"       "dataset"             "design_GSE18842"     "gene_GSE18842"       "group_GSE18842"     
dim(data_GSE18842)
#[1] 42323    91
length(unique(rownames (data_GSE18842)))
#[1] 42323 

annotation_GSE18842
#[1] "hgu133plus2.db"
str(dataset)
#chr "GSE18842"
str(design_GSE18842)
#chr "Not Paired"
str(group_GSE18842)
#'data.frame':	91 obs. of  2 variables:
#  $ Sample: Factor w/ 91 levels "GSM466947","GSM466948",..: 1 2 3 4 5 6 7 8 9 10 ...
#$ Group : chr  "d" "c" "d" "c" ...

#data  <- get(paste("gene_",dataset,sep=""))
#group <- get(paste("group_",dataset,sep=""))


##========================================================================
##  Getting the array data and annotate the probe into gene symbol
##========================================================================
dim(data_GSE18842)
colnames(data_GSE18842) == group_GSE18842$Sample 
row.names(data_GSE18842)
##  Annotate the affy probes

affyProb2Symbols <- function (aListOfProbes, db = hgu133plus2.db)
{
  require(magrittr)
  require(hgu133plus2.db)
  #require(eval(parse(text=paste0(db))))
  probeAnnoted <- AnnotationDbi::select(
    x       = db,
    keys    = aListOfProbes,
    columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
    keytype = "PROBEID"
  )
  return (probeAnnoted )
}


geneAnnoated <- affyProb2Symbols(tolower(row.names(data_GSE18842)))
#, annotation_GSE18842)
library(tidyverse)
filtered = geneAnnoated %>% distinct(PROBEID, .keep_all = TRUE)

#dim(filtered)


row.names(data_GSE18842) = filtered$SYMBOL
colnames(data_GSE18842) == group_GSE18842$Sample 


##================================
c1 <- GSA.read.gmt("X:/myGit/pathwayReview/GSA/GSA_db/c2.cp.kegg.v7.0.symbols.gmt")


clean.gmt.data <- function (db)
{
  # db <- c1 
  # str(db)
  # db$geneset.names
  # db$geneset.descriptions
  temp.geneset <- c(list())
  for (i in 1:length(db$genesets))
  {
    temp.geneset[[i]] <- db$genesets[[i]]
    if (length(which(temp.geneset[[i]]=="")) > 0 )
    {
      temp.geneset[[i]]  <- temp.geneset[[i]][-which(temp.geneset[[i]]=="")]
    }
  }
  db$genesets <- temp.geneset 
  return (db)
}

c1 <- GSA.read.gmt("X:/myGit/pathwayReview/GSA/GSA_db/Knowledge_based_gene_sets_20161005.gmt")
c1.cleaned <- clean.gmt.data (c1)
c1 <- c1.cleaned 
str(c1$genesets)


y <- rep (1, dim(group_GSE18842)[1])
y[which(group_GSE18842$Group == "d")] <- 2

gsafit <- GSA(as.matrix(data_GSE18842), y, c1$genesets, filtered$SYMBOL, resp.type = "Two class unpaired", nperms = 1000)


GSA.plot(gsafit)
GSA.listsets(gsafit, c1$geneset.names, FDRcut=.5)



str(gsafit)
plot(gsafit$GSA.scores)


str(gsafit)
plot(gsafit$GSA.scores)


gsafit$fdr.lo
gsafit$fdr.hi

# To many and will NOT print it to the html file
# gsafit$gene.scores

GSA.plot(gsafit)
GSA.listsets(gsafit, c1$geneset.names, FDRcut=.5)

GSA.genescores(5, c1$genesets, gsafit,filtered$SYMBOL)


##=======================
##  recent RNAseq data
##=======================

load ("X:/project2019/RNAseqProj/results/firstPassAnalysis/K70ALL4_DEseq2.rda")
normalizedReadCountsAll = deseq2.results.K70.all4$norm.dat

colnames(normalizedReadCountsAll )
dat.p10 <- normalizedReadCountsAll[,c(1:3, 7:8)]
y <- c(rep(1,3), rep(2,2))
gsafit.p10 <- GSA(as.matrix(dat.p10), y, c1$genesets, rownames(dat.p10), resp.type = "Two class unpaired", nperms = 1000)
GSA.plot(gsafit.p10)
GSA.listsets(gsafit.p10, c1$geneset.names, FDRcut=.5)


colnames(normalizedReadCountsAll )
dat.MD <- normalizedReadCountsAll[,c(1:6)]
y <- c(rep(1,3), rep(2,3))
gsafit.MD <- GSA(as.matrix(dat.MD), y, c1$genesets, rownames(dat.MD), resp.type = "Two class unpaired", nperms = 1000)
GSA.plot(gsafit.MD)
GSA.listsets(gsafit.MD, c1$geneset.names, FDRcut=.5)


c1$geneset.names[1]



##===========================================================
##  This is GSEA 1.2 from BroadInstitute
##===========================================================


library("utils")
library("tools")
library("dplyr")
library("GSEA")

##  simple example and it works

GSEA(input.ds = system.file('extdata', 'Leukemia_hgu95av2.gct', package = 'GSEA', mustWork = TRUE),
     input.cls = system.file('extdata', 'Leukemia.cls', package = 'GSEA', mustWork = TRUE),
     input.chip = system.file('extdata', 'Human_AFFY_HG_U95_MSigDB_7_0_final.chip',
                              package = 'GSEA', mustWork = TRUE), gs.db = system.file('extdata',
                                                                                      'h.all.v7.0.symbols.gmt', package = 'GSEA', mustWork = TRUE),
     collapse.dataset = TRUE, collapse.mode = 'max')

