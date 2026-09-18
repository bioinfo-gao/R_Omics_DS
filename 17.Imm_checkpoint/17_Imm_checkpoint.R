# 目标基因与免疫检查点基因的相关性分析 + 相关矩阵热图
#
# 输入: combined_RNAseq_counts.txt, gene.txt(免疫检查点基因列表, 单列无表头)
#       顶部 geneName 指定目标基因
# 输出: corResult.txt(显著相关的基因对), corpot.pdf(相关矩阵图)
#
# ⚠⚠ 本文件此前被一次有损编码转换严重破坏: 共 6 行代码被并进了注释,
#   结果是 group / x 未定义、不产生任何输出文件、dev.off() 关错设备。
#   六处均已还原, 每处都在原地注明了判据。

setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/17.Imm_checkpoint")
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(corrplot)

pFilter=0.05             
geneName="TSPAN6"           
expFile="combined_RNAseq_counts.txt"     
#expFile="../01.New_TCGA/combined_RNAseq_counts.txt"      
geneFile="gene.txt"       


rt=read.table(expFile, header=T, sep="\t", check.names=F)
#rt=as.matrix(rt)
rt[1:2, 1:5]

rownames(rt)=make.unique(rt[,1]) # 1
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
	
# 已还原: 原为 gene=read.table(...), 变量名被编码事故截成了 e;
# 判据是下一行使用 gene[,1], 且 27.Cuproptosis/cu.R 有逐字相同的写法。
# 取目标基因 + 免疫检查点基因的子集, 转成 样本×基因 后取 log2
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data), as.vector(gene[,1]))
data=t(data[c(geneName, sameGene),])
data=log2(data+1)

# 已还原: 一次有损编码转换把上一行的中文注释与本行代码并成了一行,
# 于是整个赋值都落进了注释里, 导致下一行的 group 未定义。
group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
# 逐个免疫检查点基因与目标基因做 Pearson 相关
# ⚠ 只有 p<pFilter 的才写入 outTab, 因此 corResult.txt 是【已筛选】的结果,
#   不能用它估计整体相关性分布, 也没有做多重检验校正。
data=t(avereps(data))

# 已还原: 原为 x=as.numeric(data[geneName,]), 同样被吞进注释,
# 导致下面 cor.test(x, y) 的 x 未定义。
x=as.numeric(data[geneName,])
outTab=data.frame()
for(i in sameGene){
	if(i==geneName){next}
    y=as.numeric(data[i,])
	corT=cor.test(x, y, method = 'pearson')
	cor=corT$estimate
	pvalue=corT$p.value
	if(pvalue<pFilter){
# 此处 data 已是 样本×(目标基因+显著相关基因), 故 cor() 得到的是基因间相关矩阵
		outTab=rbind(outTab, cbind(Query=geneName, Gene=i, cor, pvalue))
	}
}
# 已还原: 原为 write.table(...), 被吞进注释后本脚本【不产生任何输出文件】。
write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)

# 已还原: 原为 data=t(...), 把矩阵转成 样本×基因 并只保留目标基因与显著相关基因,
# 这是下一行 cor(data) 能得到基因间相关矩阵的前提。
data=t(data[c(geneName, as.vector(outTab[,2])),])
M=cor(data)

# 已还原: 原为 pdf(file="corpot.pdf",...), 被吞进注释后图不落盘,
# 且文件末尾的 dev.off() 关的是默认设备。
pdf(file="corpot.pdf",width=15,height=15)
corrplot(M,
         order="original",
         method = "color",
         number.cex = 0.7,
         addCoef.col = "black",
         diag = TRUE,
         tl.col="black",
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()


