#save PDF, dim: 9 inc x 5 inc
library("rowr")
path="/Users/GaMinh/Dropbox/WSU/Papers/PathwayReview/mmu/"
setwd(path);


Accuracy <-cbind(FE = c(0.83, 0.96, 0.95, 0.85, 0.77, 0.87, 0.78, 0.98, 0.81, 0.82, 0.95), 
                 WebGestalt = c(0.88, 0.88, 0.91, 0.73, 0.83, 0.96, 0.95, 0.86, 0.71, 0.95, 0.81),
                 GOstats = c(0.89, 0.80, 0.89, 0.82, 0.79, 0.93, 0.92, 0.87, 0.98, 0.95, NA),
                 KS = c(0.81, 0.87, 0.92, 0.85, 0.73, 0.84, 0.71, 0.98, 0.93, 0.73, 0.84), 
                 WRS = c(0.86, 0.97, 0.94, 0.85, 0.78, 0.85, 0.74, 0.97, 0.94, 0.79, 0.97), 
                 GSEA = c(0.59, 0.86, 0.84, 0.89, 0.69, 0.91, 0.75, 0.99, 0.98, 0.65, 0.98),
                 GSA = c(0.89, 0.85, 0.94, 0.83, 0.83, 0.92, 0.8, 0.88, 0.86, 0.89, 0.91), 
                 PADOG = c(0.88, 0.95, 0.95, 0.82, 0.78, 0.92, 0.86, 0.96, 0.9, 0.91, 0.92),
                 SPIA = c(0.75,0.95,0.98, 0.88, 0.79, 0.83, 0.84, 0.98, 0.73, 0.92, 0.98), 
                 ROntoTools = c(0.81, 0.97, 0.98, 0.9, 0.85, 0.91, 0.78, 0.98, 0.83, 0.94, 0.92))

pdf("Accuracy.pdf", width = 3.3, height = 6.4)
par(mar=c(6,4,4,2))
boxplot(Accuracy, las = 2, ylab = "Accuracy", notch = TRUE, 
        col=c("lightcoral","sandybrown", "darkkhaki", "olivedrab4", "green3", 
              "mediumseagreen", "chartreuse3", "mediumaquamarine", "steelblue1", "dodgerblue"))
dev.off()

Sensitivity <-cbind(FE = c(0.12, 0, 0, 0.11, 0.15, 0, 0.15, 0, 0, 0, 0), 
                    WebGestalt = c(0.00, 0.00, 0.00, 0.17, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                    GOstats = c(0.21,0.21, 0.00, 0.05, 0.00, 0.00, 0.00, 0.29, 0.00, 0.00, NA),
                    KS = c(0.65, 1, 0, 0.053, 0.21, 0, 0.2, 0, 0.5, 1, 0), 
                    WRS = c(0.71, 0, 0, 0.11, 0.21, 0, 0.25, 0, 0.5, 0, 0), 
                    GSEA = c(0.71, 0, 0.67, 0.16, 0.09, 0, 0, 0, 0, 1, 0),
                    GSA = c(0.18, 0, 0.67, 0.32, 0.18, 0, 0.1, 0, 0.5, 0, 0), 
                    PADOG = c(0.23, 0, 0.33, 0.32, 0.18, 0, 0.05, 0, 1, 0, 0), 
                    SPIA = c(0.38,1,1, 0.23, 0.33, 0, 0.1, 0, 0, 0, 0), 
                    ROntoTools = c(0.47, 1, 1, 0.26, 0.38, 0, 0.06, 0, 0, 0, 0))

pdf("Sensitivity.pdf", width = 3.3, height = 6.4)
par(mar=c(6,4,4,2))
boxplot(Sensitivity, las = 2, ylab = "Sensitivity", notch = TRUE, 
        col=c("lightcoral","sandybrown", "darkkhaki", "olivedrab4", "green3", 
              "mediumseagreen", "chartreuse3", "mediumaquamarine", "steelblue1", "dodgerblue"))
dev.off()

Specificity <-cbind(FE = c(0.9, 0.96, 0.96, 0.93, 0.89, 0.88, 0.85, 0.99, 0.81, 0.82, 0.96), 
                    WebGestalt = c(0.98, 1.00, 0.93, 0.81, 0.85, 0.97, 0.97, 0.90, 0.90, 1.00, 0.84),
                    GOstats = c(0.96, 0.93, 0.91, 0.91, 0.80, 0.94, 0.93, 0.92, 0.98, 0.97, NA),
                    KS = c(0.83, 0.87, 0.93, 0.94, 0.84, 0.85, 0.77, 0.98, 0.93, 0.73, 0.85), 
                    WRS = c(0.88, 0.97, 0.95, 0.93, 0.89, 0.86, 0.8, 0.97, 0.94, 0.8, 0.97), 
                    GSEA = c(0.58, 0.86, 0.85, 0.97, 0.82, 0.93, 0.83, 1, 0.99, 0.65, 0.99),
                    GSA = c(0.95, 0.86, 0.95, 0.89, 0.96, 0.94, 0.88, 0.89, 0.87, 0.89, 0.92), 
                    PADOG = c(0.93, 0.95, 0.95, 0.87, 0.9, 0.93, 0.94, 0.96, 0.89, 0.91, 0.93),
                    SPIA = c(0.82,0.95,0.98, 0.95, 0.91, 0.85, 0.91, 0.98, 0.73, 0.93, 0.98), 
                    ROntoTools = c(0.85, 0.97, 0.98, 0.97, 0.94, 0.92, 0.86, 0.98, 0.83, 0.94, 0.93))
pdf("Specificity.pdf", width = 3.3, height = 6.4)
par(mar=c(6,4,4,2))
boxplot(Specificity, las = 2, ylab = "Specificity", notch = TRUE, 
        col=c("lightcoral","sandybrown", "darkkhaki", "olivedrab4", "green3", 
              "mediumseagreen", "chartreuse3", "mediumaquamarine", "steelblue1", "dodgerblue"))
dev.off()

AUC <-cbind(FE = c(0.625, 0.615, 0.766, 0.495, 0.521, 0.535, 0.507, 0.76, 0.87, 0.9, 0.661),
            WebGestalt = c(0.56, 0.63, 0.61, 0.49, 0.50, 0.91, 0.50, 0.85, 0.60, 0.50, 0.66),
            GOstats = c(0.74, 0.64, 0.52, 0.55, 0.50, 0.69, 0.69, 0.87, 0.80, 0.62, 0.65),
            KS = c(0.793, 1, 0.652, 0.59, 0.635, 0.703, 0.536, 0.922, 0.616, 0.869, 0.848), 
            WRS = c(0.835, 0.940, 0.631, 0.665, 0.619, 0.588, 0.533, 0.922, 0.626, 0.785, 0.505),
            GSEA = c(0.649, 0.734, 0.867, 0.773, 0.589, 0.763, 0.652, 0.711, 0.834, 0.808, 0.765),
            GSA = c(0.792, 0.503, 0.901, 0.636, 0.785, 0.624, 0.551, 0.783, 0.775, 0.689, 0.665), 
            PADOG = c(0.814, 0.925, 0.937, 0.694, 0.623, 0.694, 0.523, 0.58, 0.923, 0.542, 0.655),
            SPIA = c(0.793,0.98,0.978, 0.719, 0.682, 0.63, 0.492, 0.755, 0.659, 0.911, 0.686),
            ROntoTools = c(0.783, 0.988, 0.981, 0.838, 0.799, 0.563, 0.634, 0.898, 0.594, 0.911, 0.794))

pdf("AUC.pdf", width = 5, height = 6.4)
par(mar=c(6,4,4,2))
boxplot(AUC, las = 2, ylab = "AUC under ROC-curve", 
        notch = TRUE, col=c("lightcoral","sandybrown", "darkkhaki", "olivedrab4", "green3", 
                            "mediumseagreen", "chartreuse3", "mediumaquamarine", "steelblue1", "dodgerblue"))
dev.off()

nonTB = as.vector(AUC[,1:8])
TB = as.vector(AUC[,9:10])

pdf("AUC_nonTBvsTB.pdf", width = 3, height = 6.4)

boxplot(nonTB, TB, names=c("Non-TB","TB"), ylab = "AUC under ROC-curve", notch = TRUE, 
        col = (c("gold","indianred3")))
# grid (NA,NULL, lty = 6, col = "cornsilk2", equilogs = FALSE) 
# boxplot(nonPT, PT, names=c("Non-TB","TB"), 
#         ylab = "AUC under ROC-curve", add = TRUE, ylim=c(0.47,1.1))
dev.off()



wilcox.test(nonTB, TB, alternative = "less")  # 0.009666
################## Classic DE genes (FC > 1.5) ###############
methods = c("FE", "KS", "WRS", "GSEA", "GSA", "PADOG", "ROntoTools", "SPIA")
# SensitivityClassic
FE = c(0,0,0,0,0,0) 
KS = c(0, 0, 0) 
WRS = c(0, 0, 0) 
GSEA = c(0.71,0,0.67,0.16,0.09,0,0,0,0,1,0)
GSA = c(0.18,0,0.67,0.32,0.18,0,0.1,0,0.5,0,0) 
PADOG = c(0.23,0,0.33,0.32,0.18,0,0.05,0,1,0,0) 
ROntoTools = c(0,0,0)
SPIA = c(1,0,0)
SensitivityClassic = cbind.fill(FE, KS, WRS, GSEA, GSA, PADOG, ROntoTools, SPIA, fill = NA)
colnames(SensitivityClassic) = methods
boxplot(SensitivityClassic)

# AccuracyClassic
FE = c(0.92, 0.99, 0.83, 0.89, 0.96, 0.99)
KS = c(0.4, 0.88, 0.99) 
WRS = c(0.4, 0.88, 0.99) 
GSEA = c(0.59,0.86,0.84,0.89,0.69,0.91,0.75,0.99,0.98,0.65,0.98)
GSA = c(0.89,0.85,0.94,0.83,0.83,0.92,0.8,0.88,0.86,0.89,0.91) 
PADOG = c(0.88,0.95,0.95,0.82,0.78,0.92,0.86,0.96,0.9,0.91,0.92)
ROntoTools = c(0.4,0.87,0.97)
SPIA = c(0.667,0.92,0.98)
AccuracyClassic = cbind.fill(FE, KS, WRS, GSEA, GSA, PADOG, ROntoTools, SPIA, fill = NA)
colnames(AccuracyClassic) = methods
boxplot(AccuracyClassic)

#SpecificityClassic
FE = c(1,1,0.99,0.99,0.96,0.99)
KS = c(1, 1, 1)
WRS = c(1, 1, 1)
GSEA = c(0.58,0.86,0.85,0.97,0.82,0.93,0.83,1,0.99,0.65,0.99)
GSA = c(0.95,0.86,0.95,0.89,0.96,0.94,0.88,0.89,0.87,0.89,0.92)
PADOG = c(0.93,0.95,0.95,0.87,0.9,0.93,0.94,0.96,0.89,0.91,0.93)
ROntoTools = c(1,0.95,0.97)
SPIA = c(0,1,0.98)
SpecificityClassic = cbind.fill(FE, KS, WRS, GSEA, GSA, PADOG, ROntoTools, SPIA, fill = NA)
colnames(SpecificityClassic) = methods
boxplot(SpecificityClassic)

#AUCClassic
FE = c(0.5, 0.998, 0.546, 0.45, 0.667, 0.582) 
KS = c(1, 0.607, 0.696)
WRS = c(1,0.621,0.696)
GSEA = c(0.649,0.734,0.867,0.773,0.589,0.763,0.652,0.711,0.834,0.808,0.765)
GSA = c(0.729,0.503,0.901,0.636,0.785,0.624,0.551,0.783,0.775,0.689,0.665)
PADOG = c(0.814,0.925,0.937,0.694,0.623,0.694,0.523,0.58,0.923,0.542,0.655)
ROntoTools = c(0.833,0.695,0.683) 
SPIA = c(1, 0.932, 0.593)
AUCClassic = cbind.fill(FE, KS, WRS, GSEA, GSA, PADOG, ROntoTools, SPIA, fill = NA)
colnames(AUCClassic) = methods
boxplot(AUCClassic)


######## ONLY DATA SETS WITH MORE THAN 2 DE GENES ############


################## Classic DE genes (FC > 1.5) ###############
methods = c("FE", "KS", "WRS", "GSEA", "GSA", "PADOG", "ROntoTools", "SPIA")
# SensitivityClassic
FE15 = c(0, 0, 0) 
KS15 = c(0, 0, 0) 
WRS15 = c(0, 0, 0) 
GSEA15 = c(0.67, 0, 0)
GSA15 = c(0.67, 0.1, 0.5) 
PADOG15 = c(0.33, 0.05, 1) 
ROntoTools15 = c(0,0,0)
SPIA15 = c(1,0,0)
SensitivityClassic15 = cbind.fill(FE15, KS15, WRS15, GSEA15, GSA15, PADOG15, ROntoTools15, SPIA15, fill = NA)
colnames(SensitivityClassic15) = methods
boxplot(SensitivityClassic15)

# AccuracyClassic
FE15 = c(0.99, 0.89, 0.96)
KS15 = c(0.4, 0.88, 0.99) 
WRS15 = c(0.4, 0.88, 0.99) 
GSEA15 = c(0.84, 0.75, 0.98)
GSA15 = c(0.94, 0.8, 0.86) 
PADOG15 = c(0.95, 0.86, 0.9)
ROntoTools15 = c(0.4, 0.87, 0.97)
SPIA15 = c(0.667,0.92,0.98)
AccuracyClassic15 = cbind.fill(FE15, KS15, WRS15, GSEA15, GSA15, PADOG15, ROntoTools15, SPIA15, fill = NA)
colnames(AccuracyClassic) = methods
boxplot(AccuracyClassic)

#SpecificityClassic
FE15 = c(1, 0.99, 0.96)
KS15 = c(1, 1, 1)
WRS15 = c(1, 1, 1)
GSEA15 = c(0.85, 0.83, 0.99)
GSA15 = c(0.95, 0.88, 0.87)
PADOG15 = c(0.95, 0.94, 0.89)
ROntoTools15 = c(1, 0.95, 0.97)
SPIA15 = c(0, 1, 0.98)
SpecificityClassic15 = cbind.fill(FE15, KS15, WRS15, GSEA15, GSA15, PADOG15, ROntoTools15, SPIA15, fill = NA)
colnames(SpecificityClassic) = methods
boxplot(SpecificityClassic)

#AUCClassic
FE15 = c(0.998, 0.45, 0.667) 
KS15 = c(1, 0.607, 0.696)
WRS15 = c(1,0.621,0.696)
GSEA15 = c(0.867, 0.652, 0.834)
GSA15 = c(0.901, 0.551, 0.775)
PADOG15 = c(0.937, 0.523, 0.923)
ROntoTools15 = c(0.833, 0.695, 0.683) 
SPIA15 = c(1, 0.932, 0.593)
AUCClassic13 = cbind.fill(FE15, KS15, WRS15, GSEA15, GSA15, PADOG15, ROntoTools15, SPIA15, fill = NA)
colnames(AUCClassic) = methods
boxplot(AUCClassic)

################## Classic DE genes (FC > 1.3) ###############
methods = c("FE", "KS", "WRS", "GSEA", "GSA", "PADOG", "ROntoTools", "SPIA")
# SensitivityClassic
FE13 = c(0, 0, 0) 
KS13 = c(0, 0, 0) 
WRS13 = c(0, 0, 0) 
GSEA13 = c(0.67, 0, 0)
GSA13 = c(0.67, 0.1, 0.5) 
PADOG13 = c(0.33, 0.05, 1) 
ROntoTools13 = c(0,0,0)
SPIA13 = c(1,0,0)
SensitivityClassic13 = cbind.fill(FE13, KS13, WRS13, GSEA13, GSA13, PADOG13, ROntoTools13, SPIA13, fill = NA)
colnames(SensitivityClassic) = methods
boxplot(SensitivityClassic)

# AccuracyClassic
FE13 = c(0.99, 0.89, 0.97)
KS13 = c(0.4, 0.89, 0.99) 
WRS13 = c(0.4, 0.88, 0.99) 
GSEA13 = c(0.84, 0.75, 0.98)
GSA13 = c(0.94, 0.8, 0.86) 
PADOG13 = c(0.95, 0.86, 0.9)
ROntoTools13 = c(0.4, 0.89, 0.97)
SPIA13 = c(0.667,0.9,0.98)
AccuracyClassic13 = cbind.fill(FE13, KS13, WRS13, GSEA13, GSA13, PADOG13, ROntoTools13, SPIA13, fill = NA)
colnames(AccuracyClassic) = methods
boxplot(AccuracyClassic)

#SpecificityClassic
FE13 = c(1, 0.99, 0.99)
KS13 = c(1, 1, 1)
WRS13 = c(1, 1, 1)
GSEA13 = c(0.85, 0.83, 0.99)
GSA13 = c(0.95, 0.88, 0.87)
PADOG13 = c(0.95, 0.94, 0.89)
ROntoTools13 = c(1, 0.96, 0.97)
SPIA13 = c(0, 1 ,0.98)
SpecificityClassic13 = cbind.fill(FE13, KS13, WRS13, GSEA13, GSA13, PADOG13, ROntoTools13, SPIA13, fill = NA)
colnames(SpecificityClassic) = methods
boxplot(SpecificityClassic)

#AUCClassic
FE13 = c(0.998, 0.534, 0.56) 
KS13 = c(1, 0.695, 0.719)
WRS13 = c(1,0.716,0.714)
GSEA13 = c(0.867, 0.652, 0.834)
GSA13 = c(0.901, 0.551, 0.775)
PADOG13 = c(0.937, 0.523, 0.923)
ROntoTools13 = c(0.833, 0.853, 0.711) 
SPIA13 = c(1, 0.54, 0.608)
AUCClassic13 = cbind.fill(FE13, KS13, WRS13, GSEA13, GSA13, PADOG13, ROntoTools13, SPIA13, fill = NA)
colnames(AUCClassic) = methods
boxplot(AUCClassic)


################## Classic DE genes (FC > 1) ###############
methods = c("FE", "KS", "WRS", "GSEA", "GSA", "PADOG", "ROntoTools", "SPIA")
# SensitivityClassic
FE10 = c(0, 0, 0) 
KS10 = c(0, 0, 0) 
WRS10 = c(0, 0, 0) 
GSEA10 = c(0.67, 0, 0)
GSA10 = c(0.67, 0.1, 0.5) 
PADOG10 = c(0.33, 0.05, 1) 
ROntoTools10 = c(0,0,0)
SPIA10 = c(1,0,0)
SensitivityClassic10 = cbind.fill(FE10, KS10, WRS10, GSEA10, GSA10, PADOG10, ROntoTools10, SPIA10, fill = NA)
colnames(SensitivityClassic10) = methods
boxplot(SensitivityClassic10)

# AccuracyClassic
FE10 = c(0.99, 0.88, 0.98)
KS10 = c(0.4, 0.89, 0.98) 
WRS10 = c(0.4, 0.88, 0.98) 
GSEA10 = c(0.84, 0.75, 0.98)
GSA10 = c(0.94, 0.8, 0.86) 
PADOG10 = c(0.95, 0.86, 0.9)
ROntoTools10 = c(0.4, 0.89, 0.97)
SPIA10 = c(0.667, 0.89, 0.99)
AccuracyClassic10 = cbind.fill(FE10, KS10, WRS10, GSEA10, GSA10, PADOG10, ROntoTools10, SPIA10, fill = NA)
colnames(AccuracyClassic10) = methods
boxplot(AccuracyClassic10)

#SpecificityClassic
FE10 = c(1, 0.98, 0.98)
KS10 = c(1, 0.99, 0.99)
WRS10 = c(1, 0.98, 0.98)
GSEA10 = c(0.85, 0.83, 0.99)
GSA10 = c(0.95, 0.88, 0.87)
PADOG10 = c(0.95, 0.94, 0.89)
ROntoTools10 = c(1, 0.96, 0.97)
SPIA10 = c(0, 0.99 ,0.99)
SpecificityClassic10 = cbind.fill(FE10, KS10, WRS10, GSEA10, GSA10, PADOG10, ROntoTools10, SPIA10, fill = NA)
colnames(SpecificityClassic10) = methods
boxplot(SpecificityClassic10)

#AUCClassic
FE10 = c(0.998, 0.544, 0.54) 
KS10 = c(1, 0.478, 0.864)
WRS10 = c(1, 0.572, 0.748)
GSEA10 = c(0.867, 0.652, 0.834)
GSA10 = c(0.901, 0.551, 0.775)
PADOG10 = c(0.937, 0.523, 0.923)
ROntoTools10 = c(0.833, 0.711, 0.59) 
SPIA10 = c(1, 0.566, 0.52)
AUCClassic10 = cbind.fill(FE10, KS10, WRS10, GSEA10, GSA10, PADOG10, ROntoTools10, SPIA10, fill = NA)
colnames(AUCClassic10) = methods
boxplot(AUCClassic10)


################## Classic DE genes (FC > 0.5) ###############
methods = c("FE", "KS", "WRS", "GSEA", "GSA", "PADOG", "ROntoTools", "SPIA")
# SensitivityClassic
FE05 = c(0, 0, 0.1, 0) 
KS05 = c(0, 0, 0.053, 0) 
WRS05 = c(0, 0, 0, 0) 
GSEA05 = c(0.67, 0, 0, 0)
GSA05 = c(0.67, 0, 0.1, 0.5) 
PADOG05 = c(0.33, 0, 0.05, 1) 
ROntoTools05 = c(0, 0, 0, 0)
SPIA05 = c(1, 0, 0, 0)
SensitivityClassic05 = cbind.fill(FE05, KS05, WRS05, GSEA05, GSA05, PADOG05, ROntoTools05, SPIA05, fill = NA)
colnames(SensitivityClassic05) = methods
boxplot(SensitivityClassic05)

# AccuracyClassic
FE05 = c(0.99, 0.98, 0.81, 0.94)
KS05 = c(0.4, 0.89, 0.83, 0.95) 
WRS05 = c(0.4, 0.89, 0.85, 0.95) 
GSEA05 = c(0.84, 0.91, 0.75, 0.98)
GSA05 = c(0.94, 0.92, 0.8, 0.86) 
PADOG05 = c(0.95, 0.92, 0.86, 0.9)
ROntoTools05 = c(0.4, 0.71, 0.83, 0.96)
SPIA05 = c(0.667, 0.99, 0.86, 0.99)
AccuracyClassic05 = cbind.fill(FE05, KS05, WRS05, GSEA05, GSA05, PADOG05, ROntoTools05, SPIA05, fill = NA)
colnames(AccuracyClassic05) = methods
boxplot(AccuracyClassic05)

#SpecificityClassic
FE05 = c(1, 0.99, 0.89, 0.95)
KS05 = c(1, 1, 0.93, 0.96)
WRS05 = c(1, 1, 0.95, 0.96)
GSEA05 = c(0.85, 0.93, 0.83, 0.99)
GSA05 = c(0.95, 0.94, 0.88, 0.87)
PADOG05 = c(0.95, 0.93, 0.94, 0.89)
ROntoTools05 = c(1, 0.83, 0.92, 0.97)
SPIA05 = c(0, 1, 0.96, 0.99)
SpecificityClassic05 = cbind.fill(FE05, KS05, WRS05, GSEA05, GSA05, PADOG05, ROntoTools05, SPIA05, fill = NA)
colnames(SpecificityClassic05) = methods
boxplot(SpecificityClassic05)

#AUCClassic
FE05 = c(0.998, 0.62, 0.508, 0.51) 
KS05 = c(1, 0.75, 0.544, 0.872)
WRS05 = c(1, 0.75, 0.591, 0.812)
GSEA05 = c(0.867, 0.763, 0.652, 0.834)
GSA05 = c(0.901, 0.624, 0.551, 0.775)
PADOG05 = c(0.937, 0.694, 0.523, 0.923)
ROntoTools05 = c(0.833, 0.833, 0.618, 0.561) 
SPIA05 = c(1, 0.5, 0.554, 0.465)
AUCClassic05 = cbind.fill(FE05, KS05, WRS05, GSEA05, GSA05, PADOG05, ROntoTools05, SPIA05, fill = NA)
colnames(AUCClassic05) = methods
boxplot(AUCClassic05)

# create a data frame
methodList = c("FE", "KS", "WRS", "GSEA", "GSA", "PADOG", "ROntoTools", "SPIA")
methods = c(rep(rep(methodList, each =3), 3), rep(methodList, each =4))
ThresholdFC = c(rep("1.5", 24), rep("1.3", 24), rep("1.0", 24), rep("0.5", 32))
Sensitivity = c(FE15, KS15, WRS15, GSEA15, GSA15, PADOG15, ROntoTools15, SPIA15,
                FE13, KS13, WRS13, GSEA13, GSA13, PADOG13, ROntoTools13, SPIA13,
                FE10, KS10, WRS10, GSEA10, GSA10, PADOG10, ROntoTools10, SPIA10,
                FE05, KS05, WRS05, GSEA05, GSA05, PADOG05, ROntoTools05, SPIA05)
data=data.frame(methods, ThresholdFC, Sensitivity)


# grouped boxplot
ggplot(data, aes(x=methods, y=Sensitivity, fill=ThresholdFC)) + 
  geom_boxplot()

# create a data frame
variety=rep(LETTERS[1:7], each=40)
treatment=rep(c("high","low"),each=20)
note=seq(1:280)+sample(1:150, 280, replace=T)
data=data.frame(variety, treatment ,  note)


# grouped boxplot
ggplot(data, aes(x=variety, y=note, fill=treatment)) + 
  geom_boxplot()





# Sensitivity
FE15 = c(0, 0, 0) 
KS15 = c(0, 0, 0) 
WRS15 = c(0, 0, 0) 
GSEA15 = c(0.67, 0, 0)
GSA15 = c(0.67, 0.1, 0.5) 
PADOG15 = c(0.33, 0.05, 1) 
ROntoTools15 = c(0,0,0)
SPIA15 = c(1,0,0)

FE13 = c(0, 0, 0) 
KS13 = c(0, 0, 0) 
WRS13 = c(0, 0, 0) 
GSEA13 = c(0.67, 0, 0)
GSA13 = c(0.67, 0.1, 0.5) 
PADOG13 = c(0.33, 0.05, 1) 
ROntoTools13 = c(0,0,0)
SPIA13 = c(1,0,0)

FE10 = c(0, 0, 0) 
KS10 = c(0, 0, 0) 
WRS10 = c(0, 0, 0) 
GSEA10 = c(0.67, 0, 0)
GSA10 = c(0.67, 0.1, 0.5) 
PADOG10 = c(0.33, 0.05, 1) 
ROntoTools10 = c(0,0,0)
SPIA10 = c(1,0,0)

FE05 = c(0, 0, 0.1, 0) 
KS05 = c(0, 0, 0.053, 0) 
WRS05 = c(0, 0, 0, 0) 
GSEA05 = c(0.67, 0, 0, 0)
GSA05 = c(0.67, 0, 0.1, 0.5) 
PADOG05 = c(0.33, 0, 0.05, 1) 
ROntoTools05 = c(0, 0, 0, 0)
SPIA05 = c(1, 0, 0, 0)

# create a data frame
methodList = c("FE", "KS", "WRS", "GSEA", "GSA", "PADOG", "ROntoTools", "SPIA")
methods = c(rep(rep(methodList, each =3), 3), rep(methodList, each =4))
ThresholdFC = c(rep("1.5", 24), rep("1.3", 24), rep("1.0", 24), rep("0.5", 32))
Sensitivity = c(FE15, KS15, WRS15, GSEA15, GSA15, PADOG15, ROntoTools15, SPIA15,
                FE13, KS13, WRS13, GSEA13, GSA13, PADOG13, ROntoTools13, SPIA13,
                FE10, KS10, WRS10, GSEA10, GSA10, PADOG10, ROntoTools10, SPIA10,
                FE05, KS05, WRS05, GSEA05, GSA05, PADOG05, ROntoTools05, SPIA05)
data=data.frame(methods, ThresholdFC, Sensitivity)


# grouped boxplot
ggplot(data, aes(x=methods, y=Sensitivity, fill=ThresholdFC)) + 
  geom_boxplot() + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) +
  scale_x_discrete(limits=methodList)



# ACCURACY
FE15 = c(0.99, 0.89, 0.96)
KS15 = c(0.4, 0.88, 0.99) 
WRS15 = c(0.4, 0.88, 0.99) 
GSEA15 = c(0.84, 0.75, 0.98)
GSA15 = c(0.94, 0.8, 0.86) 
PADOG15 = c(0.95, 0.86, 0.9)
ROntoTools15 = c(0.4, 0.87, 0.97)
SPIA15 = c(0.667,0.92,0.98)

FE13 = c(0.99, 0.89, 0.97)
KS13 = c(0.4, 0.89, 0.99) 
WRS13 = c(0.4, 0.88, 0.99) 
GSEA13 = c(0.84, 0.75, 0.98)
GSA13 = c(0.94, 0.8, 0.86) 
PADOG13 = c(0.95, 0.86, 0.9)
ROntoTools13 = c(0.4, 0.89, 0.97)
SPIA13 = c(0.667,0.9,0.98)

FE10 = c(0.99, 0.88, 0.98)
KS10 = c(0.4, 0.89, 0.98) 
WRS10 = c(0.4, 0.88, 0.98) 
GSEA10 = c(0.84, 0.75, 0.98)
GSA10 = c(0.94, 0.8, 0.86) 
PADOG10 = c(0.95, 0.86, 0.9)
ROntoTools10 = c(0.4, 0.89, 0.97)
SPIA10 = c(0.667, 0.89, 0.99)

FE05 = c(0.99, 0.98, 0.81, 0.94)
KS05 = c(0.4, 0.89, 0.83, 0.95) 
WRS05 = c(0.4, 0.89, 0.85, 0.95) 
GSEA05 = c(0.84, 0.91, 0.75, 0.98)
GSA05 = c(0.94, 0.92, 0.8, 0.86) 
PADOG05 = c(0.95, 0.92, 0.86, 0.9)
ROntoTools05 = c(0.4, 0.71, 0.83, 0.96)
SPIA05 = c(0.667, 0.99, 0.86, 0.99)

# create a data frame
methodList = c("FE", "KS", "WRS", "GSEA", "GSA", "PADOG", "ROntoTools", "SPIA")
methods = c(rep(rep(methodList, each =3), 3), rep(methodList, each =4))
ThresholdFC = c(rep("1.5", 24), rep("1.3", 24), rep("1.0", 24), rep("0.5", 32))
Accuracy = c(FE15, KS15, WRS15, GSEA15, GSA15, PADOG15, ROntoTools15, SPIA15,
             FE13, KS13, WRS13, GSEA13, GSA13, PADOG13, ROntoTools13, SPIA13,
             FE10, KS10, WRS10, GSEA10, GSA10, PADOG10, ROntoTools10, SPIA10,
             FE05, KS05, WRS05, GSEA05, GSA05, PADOG05, ROntoTools05, SPIA05)
data=data.frame(methods, ThresholdFC, Accuracy)


# grouped boxplot
ggplot(data, aes(x=methods, y=Accuracy, fill=ThresholdFC)) + 
  geom_boxplot() + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) +
  scale_x_discrete(limits=methodList)

# SPECIFICITY
FE15 = c(1, 0.99, 0.96)
KS15 = c(1, 1, 1)
WRS15 = c(1, 1, 1)
GSEA15 = c(0.85, 0.83, 0.99)
GSA15 = c(0.95, 0.88, 0.87)
PADOG15 = c(0.95, 0.94, 0.89)
ROntoTools15 = c(1, 0.95, 0.97)
SPIA15 = c(0, 1, 0.98)

FE13 = c(1, 0.99, 0.99)
KS13 = c(1, 1, 1)
WRS13 = c(1, 1, 1)
GSEA13 = c(0.85, 0.83, 0.99)
GSA13 = c(0.95, 0.88, 0.87)
PADOG13 = c(0.95, 0.94, 0.89)
ROntoTools13 = c(1, 0.96, 0.97)
SPIA13 = c(0, 1 ,0.98)

FE10 = c(1, 0.98, 0.98)
KS10 = c(1, 0.99, 0.99)
WRS10 = c(1, 0.98, 0.98)
GSEA10 = c(0.85, 0.83, 0.99)
GSA10 = c(0.95, 0.88, 0.87)
PADOG10 = c(0.95, 0.94, 0.89)
ROntoTools10 = c(1, 0.96, 0.97)
SPIA10 = c(0, 0.99 ,0.99)

FE05 = c(1, 0.99, 0.89, 0.95)
KS05 = c(1, 1, 0.93, 0.96)
WRS05 = c(1, 1, 0.95, 0.96)
GSEA05 = c(0.85, 0.93, 0.83, 0.99)
GSA05 = c(0.95, 0.94, 0.88, 0.87)
PADOG05 = c(0.95, 0.93, 0.94, 0.89)
ROntoTools05 = c(1, 0.83, 0.92, 0.97)
SPIA05 = c(0, 1, 0.96, 0.99)

# create a data frame
methodList = c("FE", "KS", "WRS", "GSEA", "GSA", "PADOG", "ROntoTools", "SPIA")
methods = c(rep(rep(methodList, each =3), 3), rep(methodList, each =4))
ThresholdFC = c(rep("1.5", 24), rep("1.3", 24), rep("1.0", 24), rep("0.5", 32))
Specificity = c(FE15, KS15, WRS15, GSEA15, GSA15, PADOG15, ROntoTools15, SPIA15,
                FE13, KS13, WRS13, GSEA13, GSA13, PADOG13, ROntoTools13, SPIA13,
                FE10, KS10, WRS10, GSEA10, GSA10, PADOG10, ROntoTools10, SPIA10,
                FE05, KS05, WRS05, GSEA05, GSA05, PADOG05, ROntoTools05, SPIA05)
data=data.frame(methods, ThresholdFC, Specificity)


# grouped boxplot
ggplot(data, aes(x=methods, y=Specificity, fill=ThresholdFC)) + 
  geom_boxplot() + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) +
  scale_x_discrete(limits=methodList)


# AUC
FE15 = c(0.998, 0.45, 0.667) 
KS15 = c(1, 0.607, 0.696)
WRS15 = c(1,0.621,0.696)
GSEA15 = c(0.867, 0.652, 0.834)
GSA15 = c(0.901, 0.551, 0.775)
PADOG15 = c(0.937, 0.523, 0.923)
ROntoTools15 = c(0.833, 0.695, 0.683) 
SPIA15 = c(1, 0.932, 0.593)

FE13 = c(0.998, 0.534, 0.56) 
KS13 = c(1, 0.695, 0.719)
WRS13 = c(1,0.716,0.714)
GSEA13 = c(0.867, 0.652, 0.834)
GSA13 = c(0.901, 0.551, 0.775)
PADOG13 = c(0.937, 0.523, 0.923)
ROntoTools13 = c(0.833, 0.853, 0.711) 
SPIA13 = c(1, 0.54, 0.608)

FE10 = c(0.998, 0.544, 0.54) 
KS10 = c(1, 0.478, 0.864)
WRS10 = c(1, 0.572, 0.748)
GSEA10 = c(0.867, 0.652, 0.834)
GSA10 = c(0.901, 0.551, 0.775)
PADOG10 = c(0.937, 0.523, 0.923)
ROntoTools10 = c(0.833, 0.711, 0.59) 
SPIA10 = c(1, 0.566, 0.52)

FE05 = c(0.998, 0.62, 0.508, 0.51) 
KS05 = c(1, 0.75, 0.544, 0.872)
WRS05 = c(1, 0.75, 0.591, 0.812)
GSEA05 = c(0.867, 0.763, 0.652, 0.834)
GSA05 = c(0.901, 0.624, 0.551, 0.775)
PADOG05 = c(0.937, 0.694, 0.523, 0.923)
ROntoTools05 = c(0.833, 0.833, 0.618, 0.561) 
SPIA05 = c(1, 0.5, 0.554, 0.465)

# create a data frame
methodList = c("FE", "KS", "WRS", "GSEA", "GSA", "PADOG", "ROntoTools", "SPIA")
methods = c(rep(rep(methodList, each =3), 3), rep(methodList, each =4))
ThresholdFC = c(rep("1.5", 24), rep("1.3", 24), rep("1.0", 24), rep("0.5", 32))
AUC = c(FE15, KS15, WRS15, GSEA15, GSA15, PADOG15, ROntoTools15, SPIA15,
        FE13, KS13, WRS13, GSEA13, GSA13, PADOG13, ROntoTools13, SPIA13,
        FE10, KS10, WRS10, GSEA10, GSA10, PADOG10, ROntoTools10, SPIA10,
        FE05, KS05, WRS05, GSEA05, GSA05, PADOG05, ROntoTools05, SPIA05)
data=data.frame(methods, ThresholdFC, AUC)


# grouped boxplot
ggplot(data, aes(x=methods, y=AUC, fill=ThresholdFC)) + 
  geom_boxplot() + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) +
  scale_x_discrete(limits=methodList)
