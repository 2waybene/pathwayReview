# library(vioplot)
library(ggplot2)
library(ggthemes)
# library(psych)
library(data.table)

############### FIGURE 2: ranks and p-values of individual method ############### 
# library(vioplot)
setwd("x:/myGit/pathwayReview/")
path="x:/myGit/pathwayReview/"

filename="x:/myGit/pathwayReview/reviewDatasets/datasetslist.txt"
rndataSets = as.character(t(read.table(filename,header=F,sep="\t",stringsAsFactors=F)))

methods <- c("FE", "WebGestalt", "GOstats", "KS", "WRS", "GSEA", "GSA", "PADOG", 
            "SPIA", "ROntoTools", "CePaGSA", "CePaORA", "PathNet") #change order
Type <- c("Non-TB", "Non-TB", "Non-TB", "Non-TB", "Non-TB", "Non-TB", "Non-TB", "Non-TB", 
          "TB", "TB", "TB", "TB", "TB")

load("WebGestaltResult.RData")
load("GOstatsResult.RData")
addition <- list(WebGestalt = WebGestaltResult, GOstats = GOstats)
load("ResultRankPValue.RData")
RankPvalue <- append(result, addition, 1)

summary <- as.data.frame(rbindlist(RankPvalue, fill=FALSE, idcol="ID"))
summary["Type"] <- rep(Type, each = 75)
colnames(summary) <- c("Method", "Rank", "Pvalue", "Type")
summary = transform(summary, Method=factor(Method,levels=methods))

wt <- lapply(RankPvalue, function(x) median(x$Rank)) #median of the method
wt <- unlist(wt)

PDFPath = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/Rank_individual.pdf"
pdf(file=PDFPath, width = 10, height = 5)

ggplot(summary, aes(x=Method, y=Rank, fill = Method)) + 
  geom_violin()+
  facet_wrap(~ Method, scale="free_x", nrow = 1, strip.position = "bottom") +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  stat_summary(fun.y=median, geom="point", size=2, color="black") +
  theme(axis.text.x = element_text(angle=45, vjust=0.6))+ 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) +
  labs(x="Methods",
       y="Ranks of target pathways") + 
  theme(legend.position = "none") + 
  geom_hline(aes(yintercept = wt[1]), data = subset(summary, Method == methods[1])) + 
  geom_hline(aes(yintercept = wt[2]), data = subset(summary, Method == methods[2])) +
  geom_hline(aes(yintercept = wt[3]), data = subset(summary, Method == methods[3])) +
  geom_hline(aes(yintercept = wt[4]), data = subset(summary, Method == methods[4])) +
  geom_hline(aes(yintercept = wt[5]), data = subset(summary, Method == methods[5])) +
  geom_hline(aes(yintercept = wt[6]), data = subset(summary, Method == methods[6])) +
  geom_hline(aes(yintercept = wt[7]), data = subset(summary, Method ==  methods[7])) +
  geom_hline(aes(yintercept = wt[8]), data = subset(summary, Method == methods[8])) +
  geom_hline(aes(yintercept = wt[9]), data = subset(summary, Method == methods[9])) +
  geom_hline(aes(yintercept = wt[10]), data = subset(summary, Method == methods[10])) +
  geom_hline(aes(yintercept = wt[11]), data = subset(summary, Method == methods[11])) +
  geom_hline(aes(yintercept = wt[12]), data = subset(summary, Method == methods[12])) +
  geom_hline(aes(yintercept = wt[13]), data = subset(summary, Method == methods[13])) 

dev.off()

wt <- lapply(RankPvalue, function(x) median(na.omit(x$pValue))) #median of the method
wt <- unlist(wt)

PDFPath = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/Pvalue_individual.pdf"
pdf(file=PDFPath, width = 10, height = 5)

ggplot(summary, aes(x=Method, y=Pvalue, fill = Method)) + 
  geom_violin()+
  facet_wrap(~ Method, scale="free_x", nrow = 1, strip.position = "bottom") +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  stat_summary(fun.y=median, geom="point", size=2, color="black") +
  theme(axis.text.x = element_text(angle=45, vjust=0.6))+ 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) +
  labs(x="Methods",
       y="P-values of target pathways") + 
  theme(legend.position = "none") + 
  geom_hline(aes(yintercept = wt[1]), data = subset(summary, Method == methods[1])) + 
  geom_hline(aes(yintercept = wt[2]), data = subset(summary, Method == methods[2])) +
  geom_hline(aes(yintercept = wt[3]), data = subset(summary, Method == methods[3])) +
  geom_hline(aes(yintercept = wt[4]), data = subset(summary, Method == methods[4])) +
  geom_hline(aes(yintercept = wt[5]), data = subset(summary, Method == methods[5])) +
  geom_hline(aes(yintercept = wt[6]), data = subset(summary, Method == methods[6])) +
  geom_hline(aes(yintercept = wt[7]), data = subset(summary, Method ==  methods[7])) +
  geom_hline(aes(yintercept = wt[8]), data = subset(summary, Method == methods[8])) +
  geom_hline(aes(yintercept = wt[9]), data = subset(summary, Method == methods[9])) +
  geom_hline(aes(yintercept = wt[10]), data = subset(summary, Method == methods[10])) +
  geom_hline(aes(yintercept = wt[11]), data = subset(summary, Method == methods[11])) +
  geom_hline(aes(yintercept = wt[12]), data = subset(summary, Method == methods[12])) +
  geom_hline(aes(yintercept = wt[13]), data = subset(summary, Method == methods[13])) 

dev.off()


############### FIGURE 3: Ranks and P-value comparison between TB and nonTB ############### 

rankNonTB <- summary[which(summary[,"Type"] == "Non-TB"),"Rank"]
rankTB <- summary[which(summary[,"Type"] == "TB"),"Rank"]



PDFPath = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/Ranks_NonTBvsTB.pdf"
pdf(file=PDFPath, width = 3, height = 5)

boxplot(rankNonTB, rankTB, col = (c("gold","indianred3")), names=c("Non-TB","TB"), 
        ylab = "Ranks of target pathways", xlab = "(a) Ranks", ylim=c(0,200), notch=TRUE)

dev.off()
wilcox.test(rankTB,rankNonTB, alternative = "less") # p-value = 0.008771
 
pvalueNonTB <- summary[which(summary[,"Type"] == "Non-TB"),"Pvalue"]
pvalueTB <- summary[which(summary[,"Type"] == "TB"),"Pvalue"]
PDFPath = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/Pvalue_NonTBvsTB.pdf"
pdf(file=PDFPath, width = 3, height = 5)

boxplot(pvalueNonTB, pvalueTB, col = (c("gold","indianred3")), names=c("Non-TB","TB"), 
        ylab = "P-values of target pathways", xlab = "(a) P-values", ylim=c(0,1.2), notch=TRUE)

dev.off()
wilcox.test(pvalueTB,pvalueNonTB, alternative = "less") # p-value = 0.000451

####### FIGURE 4 and 5: AUC, Accuracy, Specificity, Sensitivity, see OverviewMMU_KO_gene.R ####### 

######## Figure 7: Number of biased pathways calculated based on Pearson's moment coefficient

blues <- brewer.pal(9, "Blues")
blue_range <- colorRampPalette(blues)

df <- data.frame(Methods = methods, Number = b)
df = df[order(df$Number),]

PDFPath = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/NrBiasedPathway(total).pdf"
pdf(file=PDFPath, width = 11, height = 6)

ggplot(df, aes(Methods, Number)) + 
  geom_bar(stat="identity", aes(fill = Methods),colour="black")+  #for column that has zero data point
  scale_x_discrete(limits=df$Methods) + 
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = "none",
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) + 
  ylab("Number of pathways biased") +
  scale_y_continuous(limits=c(0, 150)) +
  geom_text(size = 7,aes(label=Number), position=position_dodge(width=0.9), vjust=-0.75)

dev.off()

df0 <- data.frame(Methods = methods, Number = b0)
PDFPath = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/NrBiasedPathway(bias0).pdf"
pdf(file=PDFPath, width = 5, height = 8)


ggplot(df0, aes(Methods, Number)) + 
  geom_bar(stat="identity", aes(fill = Methods),colour="black")+  #for column that has zero data point
  scale_x_discrete(limits=df$Methods) + 
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = "none",
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) + 
  ylab("Number of pathways biased toward 0") +
  scale_y_continuous(limits=c(0, 150)) +
  geom_text(size = 5,aes(label=Number), position=position_dodge(width=0.9), vjust=-0.75)

dev.off()

df1 <- data.frame(Methods = methods, Number = b1)
PDFPath = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/NrBiasedPathway(bias1).pdf"
pdf(file=PDFPath, width = 5, height = 8)


ggplot(df1, aes(Methods, Number)) + 
  geom_bar(stat="identity", aes(fill = Methods),colour="black")+  #for column that has zero data point
  scale_x_discrete(limits=df$Methods) + 
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = "none",
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) + 
  ylab("Number of pathways biased toward 1") +
  scale_y_continuous(limits=c(0, 150)) +
  geom_text(size = 5,aes(label=Number), position=position_dodge(width=0.9), vjust=-0.75)

dev.off()

####### FIGURE 8: Number of biased methods for each pathway ##########
path = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/"
setwd(path)
load("pwNames.RData")
methods <- c("FE", "WebGestalt", "GOstats", "KS", "WRS", "GSEA", "GSA", "PADOG", 
             "SPIA", "ROntoTools", "CePaGSA", "CePaORA", "PathNet") #change order

b0 = vector(mode = "numeric", length = length(methods))
b1 = vector(mode = "numeric", length = length(methods))
b = vector(mode = "numeric", length = length(methods))
names(b0) = methods
names(b1) = methods
names(b) = methods
bias0 = data.frame(Pathway = pwNames$Name, row.names = pwNames$Symbol)
bias1 = data.frame(Pathway = pwNames$Name, row.names = pwNames$Symbol)

for (i in 1:length(methods)) {
  load(paste0(methods[i], "_bias.RData"))
  b0[i] = length(bias0List)
  b1[i] = length(bias1List)
  b[i] = b0[i] + b1[i]
  x = rep(0, nrow(bias0))
  names(x) = pwNames$Name
  x[bias0List] = 1
  bias0[methods[i]] = x
  
  x = rep(0, nrow(bias1))
  names(x) = pwNames$Name
  x[bias1List] = 1
  x = x[intersect(names(x), bias1$Pathway)]
  bias1[methods[i]] = x
  
  rm(bias0List, bias1List)
}

bias0["Total"] <- apply(bias0[,-1], 1, sum)
bias0 = bias0[order(bias0$Total, decreasing = TRUE),]

bias1["Total"] <- apply(bias1[,-1], 1, sum)
bias1 = bias1[order(bias1$Total, decreasing = FALSE),]

pwOrder = rownames(bias1)
biasBar0 <- data.frame(Pathways = rownames(bias0), Number = bias0[,"Total"])
biasBar1 <- data.frame(Pathways = rownames(bias1), Number = bias1[,"Total"])

mydata = data.frame(PathwayID = c(sub("path:","",rownames(bias0)), sub("path:","",rownames(bias1))), 
                    Pathway = c(as.character(bias0$Pathway), as.character(bias1$Pathway)), 
                    BiasType = c(rep("toward 0", nrow(bias0)), rep("toward 1", nrow(bias1))),
                         Number= c(bias0$Total, bias1$Total))

myRes0 = mydata[which(mydata$BiasType == "toward 0"),"Number"]
names(myRes0) = mydata[which(mydata$BiasType == "toward 0"),"Pathway"]
myRes0 = myRes0[order(myRes0)]
myRes1 = mydata[which(mydata$BiasType == "toward 1"),"Number"]
names(myRes1) = mydata[which(mydata$BiasType == "toward 1"),"Pathway"]
myRes = as.data.frame(cbind(myRes0, myRes1[names(myRes0)]))
colnames(myRes) = c("Bias toward 0", "Bias toward 1")
myRes["Total"] <- apply(myRes, 1, sum)
myRes = myRes[order(myRes[,3], decreasing = TRUE),]

PDFPath = "/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/NrMethodsBiased.pdf"
pdf(file=PDFPath, width = 8, height = 8)

ggplot(mydata, aes(fill=BiasType, y=Number, x=Pathway)) +
  geom_bar( stat="identity") + 
  scale_x_discrete(limits=rownames(myRes))+
  scale_y_continuous(breaks=seq(0, 12, by = 2))+ 
  labs(y = "Number of methods biased", x = "Pathways")+ 
  #theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=16,face="bold"),
        legend.text=element_text(size=14), 
        axis.text.y=element_blank())+   # Centre plot title
  scale_fill_manual('Dose', values = c('deepskyblue3', 'firebrick1')) +
  #scale_colour_manual(values=c("red", "gren"))+
  #scale_fill_brewer(palette = "Dark3") +
  #scale_colour_brewer() +
  guides(fill=guide_legend(title="Bias type")) +
  coord_flip() 

dev.off()

load("/Users/minhnguyen/Dropbox/WSU/Papers/PathwayReview/KEGG65.15Pathways_Minh.RData")
myRes <- myRes[150:1,]
mappingNameIDpathway <- names(kpn)
names(mappingNameIDpathway) <- kpn
myRes["pathID"] <- mappingNameIDpathway[rownames(myRes)]
myRes["Pathway Names"] <- rownames(myRes)
rownames(myRes) <- myRes[,"pathID"]
myRes <- myRes[,c("Pathway Names", "Bias toward 0", "Bias toward 1", "Total")]
rownames(myRes) <- sub("path:","", rownames(myRes))
write.csv(myRes, file = "TableS2.csv")
