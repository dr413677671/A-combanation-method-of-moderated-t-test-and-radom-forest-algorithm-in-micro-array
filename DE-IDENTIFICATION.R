# ===========================================================================
#
#               Differentlt expressed genes annotation application
#
#  The code are written based on:
#  empritical bayes and moderated t-test
#
# ===========================================================================
#
# Author:  Ran D (413677671@qq.com)
#
# Program Features:  
#   Read and merged affy-metrix transcripts and do the quality control
#   Normalized the DE transcipts (RMA)
#   Identified DE transcripts by fitting a linear regression model (moderated t-test)
#   Retrieved the visualized results
#
# ========================

source("http://bioconductor.org/biocLite.R")
biocLite()
library(GEOquery)
library(simpleaffy)
library(GenomicRanges)
library(oligoClasses)
library(oligo)
library(oligoData)


#Read File

#The normal human data accessories are ranged from GSM2350873 to GSM2352692(60 normal; 1760 SLe patients;)
celFiles <- list.celfiles('./newdata', full.names=TRUE)
rawData <- read.celfiles(celFiles)
geneCore <- rma(rawData, target = "core")

##This two normalization methods to caculating Expression
genePS <- rma(rawData, target = "probeset")
genePS.matrix = exprs(genePS)
geneCore.matrix = exprs(geneCore)
featureData(genePS) <- getNetAffx(genePS, "probeset")
featureData(geneCore) <- getNetAffx(geneCore, "transcript")
ph = rawData@phenoData
# Transcript Better for the research


#QC
#biocLite("RColorBrewer")
library(RColorBrewer)
cols <- brewer.pal(8, "Set1") #Set color panal

#Illustated box plots from the raw data befor been normalized
#biocLite("affyPLM")
library(affyPLM)
dir.create("outcome")
setwd("./outcome")
dir.create("QC")
setwd("QC")

svg(file="boxplot_genePS.svg",width=12,height=8)
boxplot(genePS, col=cols)
dev.off()

svg(file="boxplotgeneCore.svg",width=12,height=8)
boxplot(geneCore, col=cols)
dev.off()


#Illustrate intensity graphs

svg(file="hist_genePS.svg",width=12,height=8)
hist(genePS, col=cols)
dev.off()

svg(file="hist_geneCore.svg",width=12,height=8)
hist(geneCore, col=cols)
dev.off()

#GG plot
pmexp = pm(rawData)
sampleNames = vector()
#biocLite("ggplot2")
library(ggplot2)
logs = vector()



dir.create("RawIndensity")
setwd("./RawIndensity")
#Created microarray pictures
######################CHANGE###############################################################################
for (i in 1:240)
{
name = paste("sample",i,".svg",sep="")
svg(file=name,width=12,height=8)
image(rawData[,i],main=ph@data$index[i])
dev.off()
}

# Clustering Dendrogram
setwd("..")

svg(file="Clustering_Dendrogram_for_GenePS.svg",width=12,height=8)
eset <- exprs(genePS)
distance <- dist(t(eset),method="maximum")
clusters <- hclust(distance)
plot(clusters, cex.axis=0.25)
dev.off()

svg(file="Clustering_Dendrogram_for_GeneCore.svg",width=12,height=8)
eset <- exprs(geneCore)
distance <- dist(t(eset),method="maximum")
clusters <- hclust(distance)
plot(clusters, cex.axis=0.25)
dev.off()


# Get probe indensity from CEL files
#celfiles.qc <- fitPLM(genePS)
# Visualization
#image(celfiles.qc, which=1, add.legend=TRUE)
# The chip data may contain some artificial error
#image(celfiles.qc, which=4, add.legend=TRUE)
# Box plots are also available for the package affyPLM 
# RLE (Relative Log Expression ) should be close to zero without artificial error
#RLE(getNetAffx(genePS, "probeset"), main="RLE")
# Most mideans for NUSE (Normalised Unscaled Standard Errors are supposed to be 1.
#NUSE(celfiles.qc, main="NUSE")

#Filtering features exhibiting little variation or low signal
#biocLite("pd.hta.2.0")
library(pd.hta.2.0)
# Load Main Transcript
load(paste0(path.package("pd.hta.2.0"), "/extdata/netaffxTranscript.rda"))
transcriptAnnot <- pData(netaffxTranscript)
transcriptAnnot <- transcriptAnnot[transcriptAnnot$category == 'main', ]
transcriptMain <- transcriptAnnot$transcriptclusterid
length(transcriptMain)
# 67528 Main Transcript
# Load Main Probeset
load(paste0(path.package("pd.hta.2.0"), "/extdata/netaffxProbeset.rda"))
probesetAnnot <- pData(netaffxProbeset)
probesetAnnot <- probesetAnnot[probesetAnnot$probesettype == 'main', ]
probesetMain  <- probesetAnnot$probesetid
length(probesetMain)
# 911590 Main Probes

filtered <- nsFilter(genePS, require.entrez=FALSE, remove.dupEntrez=FALSE)
filtered <- nsFilter(geneCore, require.entrez=FALSE, remove.dupEntrez=FALSE)

#Show the features filtered by the two steps(exhibiting little viriation and low signal)
filtered$filter.log$numLowVar  #Show the number of filtered features with low viriance
#
filtered$filter.log$feature.exclude #Show the number of filtered features with low signal
#
geneCore.matrix <- exprs(filtered$eset)
geneCore.matrix <- geneCore.matrix[rownames(geneCore.matrix) %in% transcriptMain,]
genePS.matrix   <- genePS.matrix[rownames(genePS.matrix) %in% probesetMain,]


#DE Genes

library(limma)
######################CHANGE###############################################################################
#ph@data[ ,2] = c("normal","normal","normal","normal","normal","SLE","SLE","SLE","SLE","SLE")
ph@data[ ,2] = rep(c("normal","SLE"),c(60,180))
colnames(ph@data)[2]="source"
ph@data
groups = ph@data$source
f = factor(groups,levels=c("normal","SLE"))
design = model.matrix(~ 0 + f)
colnames(design) = c("normal","SLE")
fit <- lmFit(geneCore.matrix, design)
contrast.matrix <- makeContrasts(SLE_nolmal = SLE - normal, levels=design)
# The estimates of the fit for the first 10 probe sets
fit$coefficients[1:10,]
normal_fits <- contrasts.fit(fit, contrast.matrix)
normal_ebFit <- eBayes(normal_fits)
normal_ebFitDiverged <- eBayes(fit)
# The names of dataset
names(normal_ebFit)
normal_ebFit$coefficients[1:10,]
# The t-statistics and p-values of the moderated t-test for the first 10 probe sets
normal_ebFit$t[1:10,]
normal_ebFit$p.value[1:10,]
# Generating a Volcano plot
setwd("..")
dir.create("DEGenes")
setwd("./DEGenes")
name = "Volcano.svg"
svg(file=name,width=12,height=8)
volcanoplot(normal_ebFit,coef=1,highlight=10)
dev.off()
# Creating probe IDs of DE genes for functional analysis
options(digits=1)

tab <- topTable(normal_ebFit,coef=1,number=1000,adjust.method="fdr",lfc=1)
tab
# P-values (last column of tab called adj.P.Val) below a threshold (in this example the threshold is set at 0.001)
topgenes=tab
topups = topgenes[topgenes[, "logFC"] >= 1, ]
topups
topdowns = topgenes[topgenes[, "logFC"] <= -1, ]
topdowns

#DE down&low
IDs.up = rownames(topups)
IDs.down = rownames(topdowns)
DEresults = decideTests(normal_ebFit,method='global',adjust.method="BH",p.value=0.05,lfc=1)
DEresultsDiverged = decideTests(normal_ebFitDiverged,method='global',adjust.method="BH",p.value=0.05,lfc=1)
DEresults
DEresultsDiverged

save.image("myfile")



topdowns

topups

DE <- topTable(normal_ebFit, number=100000, coef=1,adjust.method="BH", lfc=5)
write.table(DE,row.names=TRUE,col.names=TRUE,quote=TRUE,file="./DEIDsinfo_ifc5.txt")
write.table(DE,row.names=FALSE,col.names=FALSE,quote=FALSE,file="./DEIDs.txt")
DE <- topTable(normal_ebFit, number=100000, coef=1,adjust.method="BH", lfc=1)
write.table(DE,row.names=TRUE,col.names=TRUE,quote=TRUE,file="./DEIDsinfo_ifc1.txt")
write.table(DE,row.names=FALSE,col.names=FALSE,quote=FALSE,file="./DEIDs_ifc1.txt")
DE <- topTable(normal_ebFit, number=100000, coef=1,adjust.method="BH", lfc=2)
write.table(DE,row.names=TRUE,col.names=TRUE,quote=TRUE,file="./DEIDsinfo_ifc2.txt")
write.table(DE,row.names=FALSE,col.names=FALSE,quote=FALSE,file="./DEIDs_ifc2.txt")
DE <- topTable(normal_ebFit, number=100000, coef=1,adjust.method="BH", lfc=4)
write.table(DE,row.names=TRUE,col.names=TRUE,quote=TRUE,file="./DEIDsinfo_ifc4.txt")
write.table(DE,row.names=FALSE,col.names=FALSE,quote=FALSE,file="./DEIDs_ifc4.txt")
DE <- topTable(normal_ebFit, number=100000, coef=1,adjust.method="BH", lfc=3)
write.table(DE,row.names=TRUE,col.names=TRUE,quote=TRUE,file="./DEIDsinfo_ifc3.txt")
write.table(DE,row.names=FALSE,col.names=FALSE,quote=FALSE,file="./DEIDs_ifc3.txt")
DE_ifc5 <- topTable(normal_ebFit, number=100000, coef=1, lfc=5)
nrow(topTable(normal_ebFit, coef=1, number=10000, lfc=5))
nrow(topTable(normal_ebFit, coef=1, number=10000, lfc=4))
nrow(topTable(normal_ebFit, coef=1, number=10000, lfc=3))
nrow(topTable(normal_ebFit, coef=1, number=10000, lfc=2))
nrow(topTable(normal_ebFit, coef=1, number=10000, lfc=1))


dim(topups)
dim(topdowns)

# Heat map
geneCore.matrix.up = geneCore.matrix[(rownames(topups)),]
geneCore.matrix.down = geneCore.matrix[(rownames(topdowns)),]
geneCore.matrix.up
geneCore.matrix.down 
sampleNames = vector()
featureNames = vector()
heatlogs = vector()
###################CHANGE################################################################################################3
if (dim(topdowns)[1]!=0)
{
for (i in 1:dim(geneCore.matrix.down)[2])
{
sampleNames = c(sampleNames,rep(colnames(geneCore.matrix.down)[i],dim(topdowns)[1]))
featureNames = c(featureNames,rownames(geneCore.matrix.down[1:dim(topdowns)[1],]))
heatlogs = c(heatlogs,geneCore.matrix.down[1:dim(topdowns)[1],i])
}

heatData = data.frame(norm_logInt=heatlogs,sampleName=sampleNames,featureName=featureNames)
svg(file="Heat_Map_GeneCore_up.svg",width=12,height=8)
dataHeat = ggplot(heatData, aes(sampleName,featureName))
dataHeat + geom_tile(aes(fill=norm_logInt)) + scale_fill_gradient(low="green", high="red")
dev.off()
}

if(dim(topups)[1]!=0)
{
for (i in 1:dim(geneCore.matrix.up)[2])
{
sampleNames = c(sampleNames,rep(colnames(geneCore.matrix.up)[i],dim(topups)[1]))
featureNames = c(featureNames,rownames(geneCore.matrix.up[1:dim(topups)[1],]))
heatlogs = c(heatlogs,geneCore.matrix.up[1:dim(topups)[1],i])
}

heatData = data.frame(norm_logInt=heatlogs,sampleName=sampleNames,featureName=featureNames)
svg(file="Heat_Map_GeneCore_up.svg",width=12,height=8)
dataHeat = ggplot(heatData, aes(sampleName,featureName))
dataHeat + geom_tile(aes(fill=norm_logInt)) + scale_fill_gradient(low="green", high="red")
dev.off()
}


# Illustrate the venn plot for genes

svg(file="Venn_GeneCore_normal_SLE.svg",width=12,height=8)
vennDiagram(DEresults)
dev.off()

svg(file="Venn_GeneCore_ori.svg",width=12,height=8)
vennDiagram(DEresultsDiverged)
dev.off()

