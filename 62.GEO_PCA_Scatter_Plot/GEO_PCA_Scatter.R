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
data=t(data)
data.pca <- prcomp(data)
# 绘制主成分的碎石图
pdf(file = "screeplot.pdf",width = 10,height = 10)
screeplot(data.pca, npcs = 10, type = "lines")
dev.off()

groupcol<-ifelse(str_detect(group$group_list ,"T"), "red",
                              "blue")
groupcol=cbind(group$X,groupcol)
groupcol=as.data.frame(groupcol)
#绘制umap图
pdf(file = "umap.pdf",height = 10,width = 10)
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

