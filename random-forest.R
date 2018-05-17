# ===========================================================================
#
#               Random forest model
#
#  The code are written based on:
#  random forest algorithm
#
# ===========================================================================
#
# Author:  Ran D (413677671@qq.com)
#
# Program Features:  
#   Chose a training-set of 30% total 
#   Generated the distribution of error-rate in different tree numbers
#   Train the model (user-defined lables; features)
#   Test the model on the 30% rest of the samples
#   Retrieved the visualized results of important features
#
# ========================

library(randomForest)
data<-read.table(file="C:/Users/DR/Desktop/genecore.matrix.up11.txt")
newdata<-c()
for(i in 1:50)
{
print(data[,i])
cbind(newdata)
}
groups<-as.factor(rep(c("normal","SLE"),c(60,180)))
data<-cbind(data,groups)
data
data<-data[ , !colnames(data) %in% rownames(datadata)[1:10]]
index <- sample(2,nrow(data),replace = TRUE,prob=c(0.7,0.3))


traindata <- data[index==2,]
testdata <- data[index==1,]
set.seed(1234)
rf_ntree <- randomForest(groups~.,data=traindata,ntree=300)
plot(rf_ntree)
#ntree=25

system.time(Randommodel <- randomForest(groups ~ ., data=traindata,importance = TRUE, proximity = FALSE, ntree = 120))
print(Randommodel)
plot(Randommodel)

importance(Randommodel)
varImpPlot(Randommodel)


head(treesize(Randommodel,terminal = TRUE))
count_data <- as.data.frame(plyr::count(treesize(Randommodel, terminal = TRUE)))
head(count_data,5)
plot(count_data)
iris_pred <- predict(Randommodel,newdata=testdata)
plot(margin(Randommodel, testdata$groupss))