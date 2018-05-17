# ===========================================================================
#
#               Subgroup distinguish
#
#  The code are written based on:
#  linear classification algorithm
#
# ===========================================================================
#
# Author:  Ran D (413677671@qq.com)
#
# Program Features:  
#   Read and merged DE transcripts(can be modified to merge probeset) 
#   Sort the DE transcipts' info generated from the identification process and illustrated the plot map
#   Linear regression the data by important features
#   Retrieved the visualized classification results
#
# ========================

library(ggplot2)
data<-read.table(file="C:/Users/DR/Desktop/genecore.matrix.up11.txt")
for(i in 1:60)
{
data<-data[-1,]
}
p<-ggplot(data,aes(x=TC01000794.hg.1,y=TC02003008.hg.1))
p+geom_boxplot(col="blue",pch=16,cex=1)+geom_point(position="jitter",col=2,pch=16,cex=1)+geom_rug()
#+geom_text(aes(y = TC02003008.hg.1 + .2, label = rownames(data),cex=0.03))

positive<-vector()
negative<-vector()
for(i in 1:180)
{
if((9-data$TC01000794.hg.1[i])>0 && 9-(data$TC02003008.hg.1[i])>0)
{
negative<-c(negative,rownames(data)[i])
}
else
{
positive<-c(positive,rownames(data)[i])
}
}

write.table(negative,row.names=FALSE,col.names=FALSE,quote=FALSE,file="C:/Users/DR/Desktop/negat.txt",sep=" ")
write.table(positive,row.names=FALSE,col.names=FALSE,quote=FALSE,file="C:/Users/DR/Desktop/posit.txt",sep=" ")

groups<-vector()
for(i in rownames(data))
{
if(i %in% positive)
{
groups<-c(groups,"positive")
}
else
{
groups<-c(groups,"negative")
}
}
write.table(groups,row.names=FALSE,col.names=FALSE,quote=TRUE,file="C:/Users/DR/Desktop/grou.txt",sep=",")