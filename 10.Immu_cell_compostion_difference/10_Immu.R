# 目标基因高/低表达两组之间的免疫细胞浸润比例差异 (CIBERSORT)
#
# 输入: geneExp.txt(需含 Type 列区分 Tumor/Normal, 第 1 列为目标基因),
#       CIBERSORT-Results.txt(CIBERSORT 输出, 末 3 列为 P-value/Correlation/RMSE)
# 输出: immune.diff.pdf

# ⚠ rm(list=ls()) 会清空调用者的整个工作环境。source 本文件前请确认没有未保存的对象。
rm(list=ls())
#
#install.packages(c("vioplot", "reshape2"))
library(limma)
library(reshape2)
library(ggpubr)
library(vioplot)
library(ggExtra)

setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/10.Immu_difference")
expFile="geneExp.txt"              #
immFile="CIBERSORT-Results.txt"    #
pFilter=0.05                       #

rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
# ⚠ 硬编码取第 1 列作为目标基因; 输入文件列序一变就换了基因且无提示
gene=colnames(rt)[1]
gene

# 只保留肿瘤样本, 并把 barcode 截到患者层级以便与 CIBERSORT 结果对齐
tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

#????目??????????量????品???蟹???
data=as.data.frame(data)
data$gene=ifelse(data[,gene]>median(data[,gene]), "High", "Low")

#??ȡ????ϸ???????ļ??????????ݽ???????
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
# CIBERSORT 的 P-value 是「去卷积整体是否可信」的检验, 不是某个细胞类型的显著性;
# 这一步剔除的是去卷积本身不可靠的样本
immune=immune[immune[,"P-value"]<pFilter,]
# 去掉末 3 列(P-value / Correlation / RMSE), 只留各细胞类型的比例
immune=as.matrix(immune[,1:(ncol(immune)-3)])

#ɾ????????Ʒ
group=sapply(strsplit(row.names(immune),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
immune=immune[group==0,]
row.names(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(immune))
immune=avereps(immune)

#???ݺϲ?
sameSample=intersect(row.names(immune), row.names(data))
rt=cbind(immune[sameSample,,drop=F], data[sameSample,,drop=F])


###dim(rt)#t(rt[1:9, 1:24])

##############????????ͼ##################
#??????ת????ggplot2?????ļ?
# ⚠ 按位置丢掉倒数第 2 列(即目标基因的表达值), 只保留 High/Low 分组列。
#   依赖上面 cbind 的列顺序; 顺序一变会丢错列且不报错。
data=rt[,-(ncol(rt)-1)]
data=melt(data,id.vars=c("gene"))
colnames(data)=c("gene", "Immune", "Expression")
#????????ͼ
group=levels(factor(data$gene))
# 固定 Low 在前、High 在后, 保证箱线图与配色顺序稳定
data$gene=factor(data$gene, levels=c("Low","High"))
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="gene",
                  xlab="",
                  ylab="Fraction",
                  legend.title=gene,
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
# 每个免疫细胞类型内部做 High vs Low 比较。
# ⚠ 这里是逐细胞类型的独立检验, 未做多重检验校正。
  stat_compare_means(aes(group=gene),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#????ͼƬ
pdf(file="immune.diff.pdf", width=7, height=6)
print(boxplot)
dev.off()
