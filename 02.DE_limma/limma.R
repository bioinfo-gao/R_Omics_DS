# limma 差异表达分析 (TCGA 肿瘤 vs 癌旁) —— 本仓库 limma 路线的主脚本
#
# 输入: 01.New_TCGA 产出的 combined_RNAseq_FPKM.txt
# 输出: limmaTab.xls(全部), diffExp.xls(达阈值), diffExpLevel.xls(达阈值基因的表达),
#       rawBox.pdf / normalBox.pdf / <gene>.diff.pdf
#
# ⚠ 本文件此前有两处被编码事故破坏(diffLab 被拆成两行、一行注释被拆散),
#   以及一处单位用错(pdf width), 均已修复并在原地注明。
# 注: 同目录的 limma_by_ZG.R 是同一分析的另一份副本。

#biocLite("limma")

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma")
# install.packages(c( "gplots", "beepr", "tidyverse" ))
# install.packages("here")
#library(here) # here 包来构建跨平台的路径，它会自动处理分隔符问题。

library(data.table)
library(tidyverse)
library(ggsignif) 
library(RColorBrewer)
library(limma)
library(ggplot2)
library(ggpubr)
library(beepr)
library(gplots)
library(pheatmap)

getwd()
setwd("C:\\Users\\zhen-\\Code\\R_code\\R_For_DS_Omics\\02.差异分析//")
    

rt=read.table("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/01.New_TCGA/combined_RNAseq_FPKM.txt",sep="\t",header=T,check.names=F)

dim(rt)
head(rt)
rt[ 1:3, 1:5 ]

rt=as.matrix(rt)

rownames(rt)=rt[,1]

# ⚠ 这里从第 3 列开始取表达值。因为 01.New_TCGA 在 symbol=T 且 RNA_type=T 时
#   写出的文件带【两列注释】: 第 1 列 gene_symbol、第 2 列 RNA 类型, 第 3 列起才是样本。
#   ⚠ 而 13.edgeR 与 14.DESeq_difference 对结构相同的文件用的是 2:ncol ——
#   若那两个脚本读的确实是同一种文件, 就会把 RNA 类型列当成一个样本,
#   as.numeric 后整列变 NA, 进而使 rowMeans 过滤掉全部基因。
#   两处口径必须统一, 请按你实际文件的列结构核对。
exp=rt[ ,3:ncol(rt)] # exp=rt[ ,2:ncol(rt)]

exp[1:2, 1:5]

dimnames=list(rownames(exp), colnames(exp))

dimnames # gene names & sample names
# ⚠ 交互式帮助调用, 属调试残留
?nrow

data=matrix( as.numeric( as.matrix(exp) ),  nrow=nrow(exp),   dimnames = dimnames) # nrow=nrow(exp) <<== is NEEDED !!

data[1:3, 1:3]

# Condense a microarray data object so that values for within-array replicate probes 
# are replaced with their average.
#?avereps()
# Not needed at all here for NGS
# data=avereps(data) 

# FPKM 上的弱过滤(仅剔除全零基因)
data=data[rowMeans(data)>0,]

data2=as.data.frame(data)

data2[1:3, 1:3]
colnames(data2)

#??01A??11A???飬??????ǰ?棬?????ź???

# 按 barcode 后缀拆出肿瘤(-01A)与癌旁(-11A)两组
exp_data_T = data2 %>% dplyr::select(str_which( colnames(.),  "-01A")) #%>% can NOT be replaced !! 

nT = ncol(exp_data_T) 

nT

exp_data_N = data2 %>% dplyr::select(str_which(colnames(.), "-11A")) # %>% can NOT be replaced 



nN = ncol(exp_data_N) 

nN

# ⚠ 列顺序在此固定为「正常在前、肿瘤在后」; 下面的 class 与 Type 都按位置构造,
#   这一行顺序一旦改动, 分组标签会静默反转。
rt= cbind(exp_data_N, exp_data_T)


dim(rt)

rt[1:2,1:5]

getwd()
#write.table(rt,file = "groupout.txt",sep="\t",quote=F)
#rt=read.table("groupout.txt",sep="\t",header=T,check.names=F,row.names = 1)

#normalize????Ϊ?????ݵ?У??
# ⚠ 已修正单位: pdf() 的 width 单位是【英寸】不是像素(原注释写的是 pixes)。
#   原值 8000 相当于要生成一张约 200 米宽的 PDF。这里改为 20 英寸宽幅,
#   样本很多时可自行调大, 但请按英寸理解这个数。
pdf(file="rawBox.pdf",  width = 20)


boxplot(rt, col = "blue", xaxt = "n", outline = F) 

dev.off()
# 分位数标准化, 拉齐样本间分布



rt=normalizeBetweenArrays(rt)

pdf(file="normalBox.pdf")

#??????ͼ
boxplot(rt,col = "red",xaxt = "n",outline = F)
dev.off()

group1=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)

conNum=length(group1[group1==1])       #????????Ʒ??Ŀ
treatNum=length(group1[group1==0])     #????????Ʒ??Ŀ

#differential????????
class <- c(rep("con",conNum),rep("treat",treatNum))  
design <- model.matrix(~factor(class)+0)
colnames(design) <- c("con","treat")
# ⚠ 方向已改为 treat - con = 肿瘤 - 正常, 即 logFC > 0 表示在【肿瘤】中更高。
#   原先是 con - treat(正常 - 肿瘤), 与 13.edgeR / 14.DESeq2 / 80.Wilcox 以及
#   已统一后的 59.GEO_TCGA_Common_DEGs 方向相反, 现已统一。
#   ⚠ 若你手上已有旧版 limmaTab.xls / diffExp.xls, 其 logFC 符号与新结果相反。

# ==============================================================
df.fit <- lmFit(rt,design)
df.matrix<- makeContrasts(treat - con,levels=design)
fit<- contrasts.fit(df.fit,df.matrix)

#??Ҷ˹????
fit2 <- eBayes(fit)

#????????
allDiff=topTable(fit2,adjust='fdr',n=Inf) 

#д??????
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)


#?页?????两?????希?pvalue小??0.05??写??????
# 已还原: 变量名 diffLab 被编码事故从中间拆成两行("diff" 与 "Lab <- ..."),
# 于是 diffLab 从未被创建, 下一行 write.table(diffLab,...) 必然报 object not found。
diffLab <- allDiff[with(allDiff, ((logFC> 1 | logFC< (-1)) & adj.P.Val < 0.05 )), ]
write.table(diffLab,file="diffExp.xls",sep="\t",quote=F)

#????????????ˮƽ?????ڹ?????
# 单基因箱线图; Type 同样按位置构造, 依赖第 77 行的列顺序
diffExpLevel <- rt[rownames(diffLab),]
write.table(diffExpLevel,file="diffExpLevel.xls",sep="\t",quote=F)

#???ӻ?
gene="THBS2" 
data=t(rt[gene,,drop=F])
Type=c(rep(1,conNum), rep(2,treatNum))
exp=cbind(data, Type)
exp=as.data.frame(exp)
colnames(exp)=c("gene", "Type")
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
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
#???




# (此处原为一行中文注释, 被编码事故拆散, 已恢复为注释)
pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
print(boxplot)
dev.off()

