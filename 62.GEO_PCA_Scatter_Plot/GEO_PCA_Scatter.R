# GEO 表达矩阵的 PCA 散点图 + 碎石图
#
# 输入: geo_exp.csv 行=基因 列=样本; group.csv 需含样本名列与分组列
# 输出: screeplot.pdf, umap.pdf
#
# ⚠ 输出文件名叫 umap.pdf, 但里面画的是 PCA, 不是 UMAP。

library(readxl)
library(tidyverse)
library(GEOquery)
library(tidyverse)
library(GEOquery)
library(limma) 
library(affy)
library(stringr)

data=read.csv("geo_exp.csv")
group<-read.csv("group.csv",header=T)
rownames(group)=group[,1]
group1=group[,2:ncol(group)]
group1 = factor(group1,levels = c("T","N"))
rt=na.omit(data)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=as.data.frame(data)
# 转置成 行=样本 列=基因, 这是 prcomp 要求的方向
data=t(data)
# ⚠ prcomp 默认 center=TRUE, scale.=FALSE。对表达矩阵不做标准化时,
#   主成分会被少数高表达/高方差基因主导。若想让各基因等权, 需 scale.=TRUE。
data.pca <- prcomp(data)
# 绘制主成分的碎石图
pdf(file = "screeplot.pdf",width = 10,height = 10)
screeplot(data.pca, npcs = 10, type = "lines")
dev.off()

# ⚠ 这里引用了 group$group_list 和下一行的 group$X, 但上面读入时并未核对过列名。
#   若 group.csv 里没有这两列, 取到的是 NULL: str_detect(NULL,...) 返回 logical(0),
#   ifelse 随之返回空向量, 最终 col 为空 —— 不报错, 但所有点会画成默认颜色。
groupcol<-ifelse(str_detect(group$group_list ,"T"), "red",
                              "blue")
groupcol=cbind(group$X,groupcol)
groupcol=as.data.frame(groupcol)
#绘制umap图
pdf(file = "umap.pdf",height = 10,width = 10)
# ⚠ 颜色与标签都是「按位置」对应的: groupcol 与 group1 的行序来自 group.csv,
#   而 data.pca$x 的行序来自 geo_exp.csv 的列序。两者若不一致, 颜色和标签会
#   静默错配到别的样本上, 图看起来完全正常。
#   稳妥做法: 用 match(rownames(data.pca$x), group$样本名列) 显式对齐后再取色。
plot(data.pca$x,cex = 2.5,main = "PCA analysis", 
     col = groupcol$groupcol,
     pch =rep(16,3))
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")
# 添加标签
text(data.pca$x,labels =group1,pos = 4,offset = 0.5,cex = 0.8)
# 添加图例
legend("bottomright",title = "Sample",inset = 0.01,
       legend = c("Tumor","Normal"),
       col = c("red","blue"),
       pch = rep(16,3))
dev.off()

