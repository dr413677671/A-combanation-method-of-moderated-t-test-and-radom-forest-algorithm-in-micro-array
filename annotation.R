# ===========================================================================
#
#               Differentlt expressed genes annotation
#
#  The code are written based on:
#  HTA 2.0 annotation reference
#
# ===========================================================================
#
# Author:  Ran D (413677671@qq.com)
#
# Program Features:  
#   Read and merged DE transcripts(can be modified to merge probeset) 
#   Sort the DE transcipts' info generated from the identification process
#   Gnerated outcome by pairing the transcripts_id and the reference file
#   Retrieved the function and gene_id for each transcipt
#
# ========================

library(pd.hta.2.0)
library(stringr)
DE<-read.table(file="./DEIDsinfo_ifc2.txt")
rowname = rownames(DE)
load(paste0(path.package("pd.hta.2.0"), "/extdata/netaffxTranscript.rda"))
transcriptAnnot <- pData(netaffxTranscript)


####The IDSData can be used for David####

##Annotation

result <- vector()

for (i in rowname)
{
	result=rbind(result,transcriptAnnot[grep(pattern=i,transcriptAnnot[,1]),])
}
colnames(result) = colnames(transcriptAnnot)
write.table(result,row.names=TRUE,col.names=TRUE,quote=TRUE,file="./result.txt")

#grep(pattern="\\", x=metaChar, value=TRUE)
genes <- str_extract_all(result$geneassignment,"//[:blank:]([:digit:]*)[:blank:]//")
genes<-unlist(genes)
genes<-unique(genes)
genes<-na.omit(genes)

dele<- function(str)
{
genesno <- vector()
for(i in str)
{
if(!(is.null(i)))
{
if(i != "")
{
genesno<-c(genesno,i)}
}
}
return(genesno)
}
genes <- genes %>% lapply(function(str){str_extract_all(str,"([:digit:]*)")}) %>%unlist(.) %>%unique(.)%>% dele(.)



write.table(genes,row.names=FALSE,col.names=FALSE,quote=FALSE,file="./genes.txt")
