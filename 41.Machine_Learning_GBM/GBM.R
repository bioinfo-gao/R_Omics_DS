library(readr)
library(VIM)
library(caret)
library(rpart)
library(rpart.plot)
library(Metrics)
library(stringr)
library(rpart)
library(tibble)
library(bitops)
library(rattle)
library(rpart.plot)
library(RColorBrewer)
library(tidyverse)
library(limma)
library(pheatmap)
library(visNetwork)
library(ggpol)
library(ggplot2)
library(sparkline)
library(dplyr)
library(tidyverse)
library(caret)
library(DALEX)
library(gbm)
library(caret)
data<-read.table(file="diffGeneExp.txt",sep = "\t",header = T,check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))

set.seed(1234)
metric <- "RMSE"
myControl <- trainControl(method="cv", number=5)

# Fitting model
fitControl <- trainControl( method = "repeatedcv", number = 4, repeats = 4)
fit <- train(x=data,y=as.factor(group),  method = "gbm", trControl = fitControl,verbose = FALSE)


#绘制基因重要性梯度图
importances <- varImp(fit)
importances
importance <- as.data.frame(importances$importance)

#输入完上面那串代码后，显示的结果就是GBM结果，将他们复制出来。
#重要性筛选区域设置多少可以自己定，我这里是只要重要性不为0都可以
#删除为重要性为0的gene后重新导入
a<-importance
varimpdf <- data.frame(var = row.names(a),
                       impor = a[,1])

ggplot(varimpdf,aes(x = reorder(var,-impor), y = impor))+
  geom_col(colour = "lightblue",fill = "lightblue")+
  labs(title="Feature gene importance (Gradient Boosting Machine)", x="",y = "importance")+
  theme(plot.title = element_text(size=12,hjust=0.5))+
  theme(axis.text.x = element_text(size = 5))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(angle = 50,vjust = 0.85,hjust = 0.75))










