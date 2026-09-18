# edgeR 差异表达分析 (TCGA 肿瘤 vs 癌旁)
#
# 输入: combined_RNAseq_counts.txt (原始 counts, 来自 01.New_TCGA)
# 输出: edgerOut.xls(全部) / diffSig.xls / up.xls / down.xls
#       normalizeExp.txt / diffmRNAExp.txt / rawBox.pdf / normalBox.pdf / <gene>.diff.pdf
#
# ⚠⚠ 本文件存在一处会影响全部结果的缺陷, 见第 67-70 行与第 91 行的 FIXME。



#BiocManager::install("edgeR")
library(edgeR)
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
library(future.apply)
# ⚠ foldChange=2 在下面是直接和 logFC 比较的(diff$logFC>foldChange),
#   即实际阈值是 log2FC > 2, 也就是 4 倍变化, 不是变量名暗示的 2 倍。
foldChange=2
padj=0.05


rt=read.table("combined_RNAseq_counts.txt",sep="\t",header=T,check.names=F) #改成自己的文件名
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]
data2=as.data.frame(data)


# ⚠ 与 59/16 等脚本同一个坑: 分组靠的是这里的列拼接顺序(正常在前),
#   下面 group 向量用 rep 按位置连续填充, 不按列名匹配。顺序一改, 标签静默反转。
#以01A和11A分组，正常放前面，肿瘤放后面
exp_data_T = data2%>% dplyr::select(str_which(colnames(.), "-01A$")) # 匹配列名或用下示写法
nT = ncol(exp_data_T) 

exp_data_N = data2%>% dplyr::select(ends_with("-11A"))
nN = ncol(exp_data_N) 

data= cbind(exp_data_N, exp_data_T)


group1=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)

conNum  =length(group1[group1==1])       #正常组样品数目
treatNum=length(group1[group1==0])     #肿瘤组样品数目


group=c(rep("normal",conNum),rep("tumor",treatNum)) #按照自己的数据更改正常组和肿瘤组的数量
# ⚠ design 在本脚本中从未被使用 —— 下面走的是 exactTest 路线, 不需要设计矩阵。属残留代码。
design <- model.matrix(~group)

# ??DGEList

y <- DGEList(counts=data,group=group)#构建列表 edgeR

# filterByExpr 是 edgeR 推荐的表达量过滤; 注意第 27 行已先做过一次 rowMeans>1 过滤
keep.exprs<-filterByExpr(y,group = group)

y <- y[keep.exprs,,keep.lib.sizes=FALSE]

dim(y)

nsamples<-ncol(y)
# ⚠ brewer.pal 的 "Paired" 最多 12 色: 样本数 >12 会告警并只返回 12 色,
#   样本数 <3 会直接报错。样本量变动时这一行容易出问题。
col<-brewer.pal(nsamples,"Paired")
y <- calcNormFactors(y)#计算样本内标准化因子
y[["samples"]][["norm.factors"]]

# FIXME ⚠⚠ 以下 4 行是 edgeR/limma 手册里用来演示「标准化前后对比」的做法:
#   人为把样本 1 的 counts 乘 0.05、样本 2 的乘 5, 制造一个失衡的示例数据集。
#   它只应该用于画 rawBox/normalBox 两张对比图。
y2 <- y
y2$samples$norm.factors <- 1
y2$counts[,1] <- ceiling(y2$counts[,1]*0.05)
y2$counts[,2] <- y2$counts[,2]*5

#标准化前后箱线图
par(mfrow=c(1,2))
lcpm1<-cpm(y2,log=TRUE)

pdf(file="rawBox.pdf",width=20,height=20)
boxplot(lcpm1,las=2,col=col,main="",xaxt = "n")
title(main="A.Example:Unnormalised data",ylab="Log-cpm")
dev.off()

y2<-calcNormFactors(y2)
y2$samples$norm.factors

lcpm2<-cpm(y2,log=TRUE)

pdf(file="normalBox.pdf",width=20,height=20)
boxplot(lcpm2,las=2,col=col,main="",xaxt = "n")
title(main="B.Example:Normalised data",ylab="Log-cpm")
dev.off()

# FIXME ⚠⚠ 这里把被人为扭曲的 y2 赋回了 y, 于是第 92/94/106 行全部作用在污染数据上。
#   后果: edgerOut.xls / diffSig.xls / up.xls / down.xls / normalizeExp.txt /
#         diffmRNAExp.txt 以及最后的箱线图, 全都是在「样本1 counts×0.05、
#         样本2 counts×5」的数据上算出来的, 不是真实表达量。
#   正确写法应为 y <- estimateCommonDisp(y) —— 但本文件未作改动, 由你确认后再改。
y <- estimateCommonDisp(y)#计算普通的离散度
y <- estimateTagwiseDisp(y)#计算基因或miRNA范围内的离散度

et <- exactTest(y,pair = c("normal","tumor"))#进行精确检验

topTags(et)#输出排名靠前的差异miRNA信息

ordered_tags <- topTags(et, n=100000)#将差异信息存入列表

#剔除FDR值为NA的行

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff

# pseudo.counts 是 edgeR 在 common dispersion 下反算的等效 counts, 只适合可视化,
# 不应再喂给别的差异分析方法。
newData=y$pseudo.counts #将y中的counts信息通过管道存入newData

#将差异miRNA的logF、logCPM、PValue、FDR存入表格
write.table(diff,file="edgerOut.xls",sep="\t",quote=F)

#将差异miRNA：p值大于规定值、logFC大于规定值的logF、logCPM、PValue、FDR存入表格
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)

#将上调miRNA的logF、logCPM、PValue、FDR存入表格
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)

#将下调miRNA的logF、logCPM、PValue、FDR存入表格
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)

#将差异miRNA在tumor和normal样本中的表达值存入表格
#
normalizeExp=rbind(id=colnames(newData),newData)

write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F) #输出所有基因校正后的表达值（normalizeExp.txt）

diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)#输出差异基因校正后的表达值（diffmRNAExp.txt）


#可视化
#样品分组
gene="DCAF5"
# ⚠ 变量 data 在此被复用(原先是表达矩阵), 从这里起含义变成单基因的表达向量
data=t(newData[gene,,drop=F])
# ⚠ Type 同样按位置构造(前 conNum 个当 Normal), 依赖第 38 行的列顺序
Type=c(rep(1,conNum), rep(2,treatNum))
exp=cbind(data, Type)
exp=as.data.frame(exp)
colnames(exp)=c("gene", "Type")
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
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



