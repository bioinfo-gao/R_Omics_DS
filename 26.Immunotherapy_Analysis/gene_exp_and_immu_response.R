# 目标基因高/低表达两组的免疫治疗应答评分差异 (TCIA immunophenoscore)
#
# 输入: geneExp.txt(需含 Type 列区分 Tumor/Normal, 第 1 列为目标基因),
#       TCIA.txt(The Cancer Immunome Atlas 的 IPS 评分, 行=样本)
# 输出: 每个 IPS 指标一张 <指标名>.pdf 小提琴图
#
# 注: 同目录下的 _ZG 版本是同一分析的另一份副本。

setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/26.Immunotherapy_Analysis")
library(limma)
library(ggpubr)

tciaFile="TCIA.txt"        #（The Cancer Immunome Atlas, TCIA）data 
expFile="geneExp.txt"      #?????????ļ?

#??ȡ?????????ļ?
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
# ⚠ 硬编码取第 1 列作为目标基因
gene=colnames(rt)[1]

#ɾ????????Ʒ
# 只留肿瘤样本, 并把 barcode 截到患者层级以便与 TCIA 表对齐
tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

#????目??????????量????品???蟹???
# 以本批样本的中位数二分高低表达 —— 切点随队列变化
Type=ifelse(data[,gene]>median(data[,gene]), "High", "Low")
Type=factor(Type, levels=c("Low","High"))
data=cbind(as.data.frame(data), Type)

#??ȡTCIA?Ĵ????ļ?
ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)

#?ϲ?????
sameSample=intersect(row.names(ips), row.names(data))
ips=ips[sameSample, , drop=F]
data=data[sameSample, "Type", drop=F]
data=cbind(ips, data)

#???ñȽ???
# ⚠ 这一行算出的 group 会在两行之后被重新算一遍(第 35 行), 属冗余
group=levels(factor(data$Type))
data$Type=factor(data$Type, levels=c("Low", "High"))
group=levels(factor(data$Type))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#???ӻ?
# 遍历除末列 Type 之外的所有 IPS 指标
# ⚠ 每个指标独立做一次检验, 未做多重检验校正
for(i in colnames(data)[1:(ncol(data)-1)]){
	rt=data[,c(i, "Type")]
	gg1=ggviolin(rt, x="Type", y=i, fill = "Type", 
	         xlab="", ylab=i,
	         legend.title=gene,
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons)
	         #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	
	pdf(file=paste0(i, ".pdf"), width=4.8, height=4.25)
	print(gg1)
	dev.off()
}



