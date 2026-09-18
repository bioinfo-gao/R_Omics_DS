# GBM (梯度提升树) 特征基因重要性
#
# 输入: diffGeneExp.txt 行=基因, 列=样本(样本名形如 <前缀>_<组别>)
# 输出: 无落盘文件 —— 最后的 ggplot 只打印到默认设备, 未 ggsave。
#       第 45 行原注释也说明了这是「把结果复制出来」的手工流程。

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
# ⚠ 同 40: 用下划线取组别; 不匹配时会静默退化成「每样本一类」
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))

set.seed(1234)
# ⚠ metric 与下一行的 myControl 在本文件中从未被使用(下面用的是 fitControl),
#   且 RMSE 是回归指标 —— 这两行是从回归模板抄来的残留, 与实际的分类任务无关。
metric <- "RMSE"
myControl <- trainControl(method="cv", number=5)

# Fitting model
fitControl <- trainControl( method = "repeatedcv", number = 4, repeats = 4)
# y 传的是 factor, 故 train() 按分类处理, 选优指标是 Accuracy 而非上面写的 RMSE
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










