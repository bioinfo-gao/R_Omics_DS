# DESeq2 差异表达分析 (TCGA 肿瘤 vs 癌旁)
#
# 输入: ../13.edgeR/combined_RNAseq_counts.txt (注意读的是 13.edgeR 目录下的副本)
# 输出: DESeq2.diff.tsv / group.txt / rawBox.pdf / normalBox.pdf /
#       histold.pdf / histnew.pdf / <gene>.diff.pdf
#
# ⚠ logFC 方向: 下面 relevel(ref="Normal") 使 logFC>0 表示在【肿瘤】中更高,
#   这与 59.GEO_TCGA_Common_DEGs 里 TCGA 部分(con-treat)的方向相反, 合并结果时须换算。

setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/14.DESeq_difference")
#BiocManager::install("DESeq2") # there will be questions in popup, choose NO!
#BiocManager::install("DESeq2", version = "1.48.1")
library("DESeq2")
library(tidyverse)
library(ggsignif) 
library(RColorBrewer)
library(limma)
library(ggplot2)
library(ggpubr)
library(beepr)
library(gplots)
library(pheatmap)


# 读取表达量的表格
rt <- read.table( "../13.edgeR/combined_RNAseq_counts.txt",header=T,sep="\t",comment.char="",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]
data2=as.data.frame(data)


exp_data_T = data2%>% dplyr::select(str_which(colnames(.), "-01A$")) # 匹配列名或用下示写法
nT = ncol(exp_data_T) 
exp_data_N = data2%>% dplyr::select(ends_with("-11A"))
nN = ncol(exp_data_N) 
# ⚠ 列顺序固定为「正常在前, 肿瘤在后」; 下面 Type 用 rep 按位置构造,
#   不按列名匹配。这一行顺序一改, 分组标签会静默反转。
data= cbind(exp_data_N, exp_data_T)


group1=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)
conNum=length(group1[group1==1])       #正常组样品数
treatNum=length(group1[group1==0])     #肿瘤组样品数

# avereps 对重复基因取均值后会产生小数, DESeq2 要求整数 counts, 故此处向下取整
count <- floor(data)#下取
#count <- ceiling(count)#上取
# 预处理，过滤低丰度的数据
countData <- count[apply(count, 1, sum) > 0 , ]
# 读取样本分组信息
data=colnames(countData)
Type=c(rep(1,conNum), rep(2,treatNum))
exp=cbind(data, Type)
exp=as.data.frame(exp)
colnames(exp)=c("id", "Type")
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
exp=as.matrix(exp)
rownames(exp)=exp[,1]
exp=exp[,2:ncol(exp)]
exp=as.data.frame(exp)
colnames(exp)=c("condition")
write.table(exp, file="group.txt",sep="\t",quote=F)
colData <- read.table( "group.txt",header=T,sep="\t",row.names=1,comment.char="",check.names=F)

# 构建DESeq2中的对象
dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ condition)

# 指定哪一组作为对照组
dds$condition <- relevel(dds$condition, ref = "Normal")

#计算每个样本的归一化系数
# ⚠ 第 67/70/73 行这三步都会被下面第 74 行的 DESeq() 重新做一遍
#   (DESeq 内部依次调用 estimateSizeFactors / estimateDispersions / nbinomWaldTest),
#   因此它们是冗余的, 且 DESeq() 会提示 "replacing previous dispersion estimates"。
dds <- estimateSizeFactors(dds)

#估计基因的离散度
dds <- estimateDispersions(dds)

#差异分析
dds <- nbinomWaldTest(dds)
dds <- DESeq(dds)
# results() 未指定 contrast, 取设计里最后一个变量的非参照水平 vs 参照 = Tumor vs Normal
res <- results(dds)
write.table(res,"DESeq2.diff.tsv",sep="\t",quote=F,col.names = NA)

#标准化(减小离散度)
# vst 已是方差稳定后的近似 log2 尺度
vsd <- vst(dds, blind = FALSE)  ## 使用vst函数进行标准化
countData_old<- assay(dds)
countData_new<- assay(vsd)
n.sample=ncol(countData)#样本数
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))

#标准化前图
pdf(file="rawBox.pdf")
boxplot(countData_old, col = cols,main="expression value",xaxt = "n")
dev.off()
#标准化后图
pdf(file="normalBox.pdf")
boxplot(countData_new, col = cols,main="expression value",xaxt = "n")
dev.off()

pdf(file="histold.pdf")
hist(countData_old)
dev.off()

pdf(file="histnew.pdf")
hist(countData_new)
dev.off()

#可视化
gene="DCAF5"
data=t(countData_new[gene,,drop=F])
Type=c(rep(1,conNum), rep(2,treatNum))
exp=cbind(data, Type)
exp=as.data.frame(exp)
colnames(exp)=c("gene", "Type")
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
# ⚠ countData_new 来自 vst(), 已经是近似 log2 尺度; 这里又做一次 log2(x+1),
#   相当于对同一数据取了两次对数。箱线图的 y 轴不再是任何标准的表达量刻度,
#   组间差异会被压缩。若想画 VST 值请直接用, 若想画原始尺度请用 counts。
exp$gene=log2(exp$gene+1)
group=levels(factor(exp$Type))
exp$Type=factor(exp$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
boxplot=ggboxplot(exp, x="Type", y="gene", color="Type",
                  xlab="",
                  ylab=paste0(gene, " expression"),
                  legend.title="Type",
                  palette = c("blue","red"),
                  add = "jitter")+ 
  stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
print(boxplot)
dev.off()
