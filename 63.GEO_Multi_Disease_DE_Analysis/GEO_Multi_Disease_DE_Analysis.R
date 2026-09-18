# GEO 多疾病数据集的下载/注释/差异分析 (GSE75380)
#
# 输入: 联网下载 GSE75380 与平台文件 GPL13497.soft.gz
# 输出: geo_exp.csv, deg_all.txt, upanddown.csv
#
# ⚠⚠ 第 49-52 行的设计矩阵与 coef 不匹配, 会使差异分析结果失去意义, 见该处说明。
# 注: 本文件与 56.GEO_Database_Processing 是同一模板, 但分组方式、是否取 log2、
#     以及设计矩阵写法都不一样, 两者结果不可直接比较。

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
# ⚠ 与 56 相同的问题: match 对表达矩阵中不存在的探针返回 NA, exp[NA,] 产生整行 NA,
#   要到第 33 行 na.omit 才被清掉。
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
# ⚠ 分组取 title 的第一个空格前单词, 完全依赖该数据集的标题写法
group_list=str_split(pdata$title,' ',simplify = T)[,1]
table(group_list)
# ⚠ 未指定 levels, 因此分组顺序按【字母序】确定 —— 这直接决定下面 coef 指向谁
group_list <- factor(group_list,ordered = F)
table(group_list)
#####进行数据矫正#####
# ⚠ 与 56 不同: 这里【没有】取 log2。若原矩阵不是 log 尺度, 后续 logFC 的含义会不同。
#   请用上下两行的 range() 确认尺度。
gene_exp_matrix_noemal=normalizeBetweenArrays(gene_exp_matrix)
range(gene_exp_matrix_noemal)
gene_exp_matrix_noemal=na.omit(gene_exp_matrix_noemal)#去除NA列
write.csv(gene_exp_matrix_noemal,file = "geo_exp.csv")
range(gene_exp_matrix_noemal)

#####进行差异分析#####
# ⚠⚠ 这里用的是【无截距】设计 ~0+group_list: 每个系数代表的是【该组的平均表达量】,
#   而不是组间差异。因此下一步 topTable(coef=2) 检验的是「第 2 组均值是否等于 0」——
#   对表达量而言这几乎总是极显著, 得到的 deg 表不是差异表达结果。
#   多组比较的正确写法是配合 makeContrasts, 例如:
#     cont <- makeContrasts(group_listB - group_listA, levels=design)
#     fit2 <- eBayes(contrasts.fit(fit, cont)); topTable(fit2, ...)
#   或者改用带截距的 ~group_list(即 56 的写法), 此时 coef=2 才是一个真正的对比。
#   未自动修改: 本数据集有多个疾病分组, 究竟要比哪两组必须由你指定。
design=model.matrix(~0+group_list)
fit=lmFit(gene_exp_matrix_noemal,design)
fit=eBayes(fit)
# ⚠ 承上: 在无截距设计下这不是一个组间对比
deg=topTable(fit,coef=2,number = Inf)
write.table(deg, file = "deg_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
logFC=1
adj.P.Val = 0.05
k1 = (deg$adj.P.Val < adj.P.Val)&(deg$logFC < -logFC)
k2 = (deg$adj.P.Val < adj.P.Val)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
write.csv(deg,file="upanddown.csv")
