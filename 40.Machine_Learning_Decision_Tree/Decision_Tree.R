install.packages("VIM")
install.packages("Metrics")
install.packages("ggpol")
install.packages("visNetwork")
install.packages("sparkline")
#本代码未设置验证集
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


data<-read.table(file="diffGeneExp.txt",sep = "\t",header = T,check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
dim(data)
colnames(data)
aggr(data)
set.seed(123)
data2<-as.data.frame(data)
#建模
mod1<-rpart(as.factor(group)~.,data = data2,method = "class",cp=0.000001)
#显示重要性
importances <- varImp(mod1)
importances %>%
  arrange(desc(Overall))

mod1$cp
#查看模型最低CP值
#Complexity parameter是决策树每一次分裂时候最小的提升量，用于平衡模型精确度于复杂度
plotcp(mod1)
#模型优化（取最低CP值）
mod1<-rpart(as.factor(group)~.,data = data2,method = "class",cp=0.00028)
visTree(mod1,main = "Decision Tree",height = "600px",
        colorY = c("greenYellow","hotPink","yellow"),legendWidth=0.2,legendNcol=2,  nodesFontSize = 16,edgesFontSize = 10,)

