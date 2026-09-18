###加载R包
library(readxl)
library(tidyverse)
library(GEOquery)
library(tidyverse)
library(GEOquery)
library(limma) 
library(affy)
library(stringr)
###下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE75380', destdir=".", AnnotGPL = T, getGPL = T)
class(gset)
gset[[1]]

#读取平台文件
GPL_data<- getGEO(filename ="GPL13497.soft.gz", AnnotGPL = T)
GPL_data_11 <- Table(GPL_data)

#提取表达量
exp <- exprs(gset[[1]])
probe_name<-rownames(exp)

#转换ID
loc<-match(GPL_data_11[,1],probe_name)
probe_exp<-exp[loc,]
raw_geneid<-(as.matrix(GPL_data_11[,"GENE_SYMBOL"]))
index<-which(!is.na(raw_geneid))
geneid<-raw_geneid[index]
exp_matrix<-probe_exp[index,]
geneidfactor<-factor(geneid)
gene_exp_matrix<-apply(exp_matrix,2,function(x) tapply(x,geneidfactor,mean))
rownames(gene_exp_matrix)<-levels(geneidfactor)
gene_exp_matrix=na.omit(gene_exp_matrix)

#####读取分组信息#####
pdata <- pData(gset[[1]])
group_list=str_split(pdata$title,' ',simplify = T)[,1]
table(group_list)
group_list <- factor(group_list,ordered = F)
table(group_list)
#####进行数据矫正#####
gene_exp_matrix_noemal=normalizeBetweenArrays(gene_exp_matrix)
range(gene_exp_matrix_noemal)
gene_exp_matrix_noemal=na.omit(gene_exp_matrix_noemal)#去除NA列
write.csv(gene_exp_matrix_noemal,file = "geo_exp.csv")
range(gene_exp_matrix_noemal)

#####进行差异分析#####
design=model.matrix(~0+group_list)
fit=lmFit(gene_exp_matrix_noemal,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
write.table(deg, file = "deg_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
logFC=1
adj.P.Val = 0.05
k1 = (deg$adj.P.Val < adj.P.Val)&(deg$logFC < -logFC)
k2 = (deg$adj.P.Val < adj.P.Val)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
write.csv(deg,file="upanddown.csv")
