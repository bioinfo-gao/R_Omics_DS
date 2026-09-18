# XGBoost 特征基因重要性
#
# 输入: diffGeneExp.txt 行=基因, 列=样本(样本名形如 <前缀>_<组别>)
# 输出: 无落盘文件 —— 第 39 行的 plot 与末尾的 ggplot 都只打印到默认设备。
#       要留存需 ggsave(...) 或 pdf()/dev.off() 包起来。
#
# 注: 本文件与 41.Machine_Learning_GBM 是同一模板的两个副本, 只换了 method。

# ⚠ 下面这段 library 列表是从 41.GBM 整段复制来的: tidyverse/caret 各出现 3 次,
#   rpart、rpart.plot 各 2 次, 且 rpart/rattle/visNetwork/ggpol/sparkline 等
#   在本文件中根本用不到。加载顺序还会影响 dplyr 与 MASS 之间的函数遮蔽。
library(xgboost)
library(caret)
library(tidyverse)
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

set.seed(123)
data<-read.table(file="diffGeneExp.txt",sep = "\t",header = T,check.names=F, row.names=1)
data=t(data)
# ⚠ 用下划线取组别; 若样本名不含 "_" 会静默退化成「每样本一类」, 先 table(group) 确认
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
# Fitting model(用caret实现)
TrainControl <- trainControl( method = "repeatedcv", number = 10, repeats = 4)
# y 传 factor, 故按分类训练; 10 折 ×4 次重复交叉验证, 由顶部 set.seed(123) 固定
model<- train(x=data,y=as.factor(group),  method = "xgbTree", trControl = TrainControl,verbose = FALSE)


# ⚠ varImp(model) 在这里被算了两次(第 39、40 行), 对大模型是明显的重复开销
plot(varImp(model))
importance <- varImp(model)
head(importance)
important <- as.data.frame(importance$importance) 
a<-important
varimpdf <- data.frame(var = row.names(a),
                       impor = a[,1])


ggplot(varimpdf,aes(x = reorder(var,-impor), y = impor))+
  geom_col(colour = "lightblue",fill = "lightblue")+
  labs(title="Feature gene importance (XGBoost)", x="",y = "importance")+
  theme(plot.title = element_text(size=12,hjust=0.5))+
  theme(axis.text.x = element_text(size = 3))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(angle = 50,vjust = 0.85,hjust = 0.75))



