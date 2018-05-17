# ===========================================================================
#
#               GO/KEGG annotation tool
#
#  Introduction:
#  A tool for GO/KEGG annotaion based on the DE genes from previous analysis.
#
# ===========================================================================
#
# Author:  Ran D (413677671@qq.com)
#
# Program Features:  
#   Read and merged DE transcripts(can be modified to merge probeset) 
#   Sort the DE transcipts
#   GO/KEGG annotation
#   Visualized the result by plot map, network plot and flow chart.
#
# ========================

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("clusterProfiler")
library(clusterProfiler)
biocLite("topGO")
library(topGO)
biocLite("Rgraphviz")
library(Rgraphviz)
data<-read.table(file="C:/Users/DR/Desktop/genes_lfc1.txt")
data

#GO Annotation
gene<-data
ego=enrichGO(OrgDb="org.Hs.eg.db", gene = gene$V1,pvalueCutoff = 0.01,readable=TRUE)
write.csv(summary(ego),"G-enrich_lFC1.csv",row.names =F)

#Visualization
svg(file="Visualization.svg",width=12,height=8)
dotplot(ego, showCategory=30)
dev.off()
svg(file="EnrichMAP.svg",width=12,height=8)
enrichMap(ego, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
dev.off()
svg(file="Flowchart.svg",width=12,height=8)
plotGOgraph(ego)
dev.off()
svg(file="Barplot.svg",width=12,height=8)
barplot(ego,drop=TRUE,showCategory = 12)
dev.off()

#Gene Set Enrichment Analysis
gsecc <- gseGO(geneList=gene$V1, ont="CC", OrgDb=org.Hs.eg.db, verbose=F)
head(summary(gsecc))
write.csv(summary(gsecc),"GeneSetEnrich.csv",row.names =F)

#KEGG Annot
kk<-enrichKEGG(gene=gene$V1,organism='human',pvalueCutoff=0.05)
write.csv(as.matrix(kk@result),"KEGG_lFC1.csv",row.names =F)
svg(file="KEGG_Barplot.svg",width=12,height=8)
barplot(kk, drop=TRUE, showCategory=30)
dev.off()
svg(file="KEGG_Visualization.svg",width=12,height=8)
dotplot(kk, showCategory=30)
dev.off()
svg(file="KEGG_EnrichMAP.svg",width=12,height=8)
enrichMap(kk, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
dev.off()
svg(file="KEGG_Flowchart.svg",width=12,height=8)
plotKEGGgraph(kk)
dev.off()



#Compare Cluster
compare <- compareCluster(gene$V1, fun="enrichKEGG",organism="hsa", pvalueCutoff=0.05)
svg(file="CompareCluster.svg",width=12,height=8)
Plot(Compare)
dev.off()


#KEGG Pathway Vidualization
require(pathview)

Hsa04110 <- pathview( gene.data = geneList, pathway.id = ¡±hsa04110¡±, species=¡±hsa¡±, limit=list(gene = max(abs(geneList))),cpd=1)