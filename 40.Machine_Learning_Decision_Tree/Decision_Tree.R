# 决策树 (rpart) 特征基因建模与可视化
#
# 输入: diffGeneExp.txt 行=基因, 列=样本(样本名形如 <前缀>_<组别>)
# 输出: 无落盘文件 —— visTree 产生的是交互式 HTML widget, 只在 Viewer 里显示。
#       要留存需 visSave(widget, "tree.html") 或改用 rpart.plot 输出 pdf。

# 依赖安装已注释: 直接 source 本文件不再触发联网安装与 CRAN 镜像选择。缺包时手动执行。
# install.packages("VIM")
# install.packages("Metrics")
# install.packages("ggpol")
# install.packages("visNetwork")
# install.packages("sparkline")
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
# ⚠ 这里用【下划线】取组别(\\2), 而 28/29/30 用的是【短横线第 5 段】—— 仓库内约定不统一。
#   若样本名不含 "_", gsub 不匹配会原样返回整个样本名, 每个样本自成一类。
#   先 table(group) 确认类别数再往下跑。
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
dim(data)
colnames(data)
aggr(data)
set.seed(123)
data2<-as.data.frame(data)
#建模
# cp 先设得极小, 让树长到最大, 再由下面的 plotcp/cptable 挑剪枝点
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
# cp 由上一步的 cptable 自动取「交叉验证误差 xerror 最小」处, 不再手抄常数。
# (原为硬编码 cp=0.00028, 换数据集即失效且不报错)
bestcp <- mod1$cptable[which.min(mod1$cptable[,"xerror"]), "CP"]
mod1<-rpart(as.factor(group)~.,data = data2,method = "class",cp=bestcp)
# ⚠ 全部数据既用于建树也用于看效果, 文件开头也注明了「本代码未设置验证集」,
#   因此这棵树的表现不能当作泛化能力的证据。
visTree(mod1,main = "Decision Tree",height = "600px",
        colorY = c("greenYellow","hotPink","yellow"),legendWidth=0.2,legendNcol=2,  nodesFontSize = 16,edgesFontSize = 10,)

