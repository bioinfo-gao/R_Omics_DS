# 从 GEO 下载并整理表达矩阵, 做探针到基因的注释, 再跑 limma 差异分析
#
# 输入: 联网下载 GSE205185 与平台文件 GPL21185.soft.gz
# 输出: geo_exp.csv(供 59/62/63 使用), deg_all.txt, upanddown.csv
#
# ⚠ logFC 的方向见第 40 行说明, 与 59.GEO_TCGA_Common_DEGs 的约定相反。

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
gset = getGEO('GSE205185', destdir=".", AnnotGPL = T, getGPL = T)
class(gset)
gset[[1]]

#读取平台文件
GPL_data<- getGEO(filename ="GPL21185.soft.gz", AnnotGPL = T)
GPL_data_11 <- Table(GPL_data)

#提取表达量
exp <- exprs(gset[[1]])
probe_name<-rownames(exp)

#转换ID
# ⚠ match 对平台文件里存在、但表达矩阵中没有的探针返回 NA,
#   exp[NA,] 会生成整行 NA。下面第 27 行是按【基因符号是否为 NA】过滤,
#   并不能去掉这些行 —— 它们要到第 33 行的 na.omit 才被清掉。
#   更稳妥: loc 先剔除 NA, 或改用 intersect 对齐。
loc<-match(GPL_data_11[,1],probe_name)
probe_exp<-exp[loc,]
raw_geneid<-(as.matrix(GPL_data_11[,"GENE_SYMBOL"]))
index<-which(!is.na(raw_geneid))
geneid<-raw_geneid[index]
exp_matrix<-probe_exp[index,]
geneidfactor<-factor(geneid)
# 同一基因的多个探针取均值合并; 也有取中位数或最大值的做法, 口径不同结果会略有差异
gene_exp_matrix<-apply(exp_matrix,2,function(x) tapply(x,geneidfactor,mean))
rownames(gene_exp_matrix)<-levels(geneidfactor)
gene_exp_matrix=na.omit(gene_exp_matrix)

#####读取分组信息#####
pdata <- pData(gset[[1]])
# ⚠ 分组规则依赖该数据集 source_name_ch1 里的具体措辞, 换数据集必须改这一行
group_list <- ifelse(str_detect(pdata$source_name_ch1,"primary breast tumour"), "T",
                     "N")
group_list
# ⚠ 方向: levels=c("T","N") 令【T 为参照】, 下面 topTable(coef=2) 取的是 N vs T,
#   因此本脚本 logFC > 0 表示在【正常】中更高。
#   而 59.GEO_TCGA_Common_DEGs 的 GEO 部分用 levels=c("N","T"), 方向正好相反。
#   两个脚本的 deg 表不能直接混用, 合并前需换算符号。
group_list = factor(group_list,
                    levels = c("T","N"))
group_list
pdata$group=group_list

#####进行数据矫正#####

# ⚠ 这两张标准化前后的箱线图【没有落盘】: 没有 pdf() 包裹,
#   下面的 dev.off() 关掉的是隐式打开的默认设备。要留存需自行加 pdf()/dev.off()。
boxplot(gene_exp_matrix,outline=T, notch=T,col=group_list, las=2)
dev.off()
gene_exp_matrix_noemal=normalizeBetweenArrays(gene_exp_matrix)
boxplot(gene_exp_matrix_noemal,outline=T, notch=T,col=group_list, las=2)
range(gene_exp_matrix_noemal)
# ⚠ 先 normalizeBetweenArrays 再取 log2: Agilent/GEO 的表达矩阵很多已经是 log2 尺度,
#   若如此则这里等于取了两次对数。上下两行的 range() 正是用来核对尺度的 —— 请务必看一眼:
#   若标准化后的 range 已在 0~20 量级, 说明原本就是 log 尺度, 不应再 log2。
gene_exp_matrix_noemal <- log2(gene_exp_matrix_noemal+1)
gene_exp_matrix_noemal <-as.data.frame(gene_exp_matrix_noemal)
gene_exp_matrix_noemal=na.omit(gene_exp_matrix_noemal)#去除NA列
write.csv(gene_exp_matrix_noemal,file = "geo_exp.csv")
range(gene_exp_matrix_noemal)
dev.off()

#####进行差异分析#####
design=model.matrix(~group_list)
# 原注释提醒: 分组向量的顺序必须与表达矩阵的列顺序一致 —— 这里是按位置对应的, 不做名字校验
fit=lmFit(gene_exp_matrix_noemal,design)#这里要注意，分组的样本与矩阵样本是否相符，不相符则去文件中调整然后读入
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
