# Wilcoxon 秩和检验做差异表达 (肿瘤 vs 癌旁), 不依赖分布假设
#
# 输入: combined_RNAseq_TPM.txt (行=基因, 列=样本)
# 输出: wilcoxout.tsv —— |log2FC|>1 且 FDR<0.05 的基因
#
# ⚠ 方法学提醒见第 9 行: 输入是 TPM, 但下面用的是为 counts 设计的 edgeR 流程。

#读取TPM矩阵
readCount<-read.table(file="combined_RNAseq_TPM.txt", header = T, row.names = 1, stringsAsFactors = F,check.names = F)
# edger 标准化并删除低表达基因
library(edgeR)
library(limma)
# 从 barcode 第 4 段首字符取样本类型: 0=肿瘤, 1=正常(2=复发并入 1)
group1=sapply(strsplit(colnames(readCount),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)
# ⚠ 这里把 TPM 当作 counts 传给 DGEList。filterByExpr 与 calcNormFactors(TMM)
#   都是按【原始计数】的统计性质设计的: TPM 已经做过长度与深度归一化,
#   再做一次 TMM 属于重复校正, 过滤阈值的含义也不再成立。
#   若要走这条流程, 应改用 combined_RNAseq_counts.txt;
#   若坚持用 TPM, 则应跳过 DGEList/TMM, 直接对 TPM 做 Wilcoxon。
#   未自动修改: 换输入文件会改变全部结果, 需你决定。
y <- DGEList(counts=readCount,group=group1)
#删除过低表达量基因
keep <- filterByExpr(y)
y <- y[keep,keep.lib.sizes=FALSE]
##进行TMM标准化并转移到CPM（百万计数）
y <- calcNormFactors(y,method="TMM")
count_norm=cpm(y)
count_norm<-as.data.frame(count_norm)

# 对每个基因进行Wilcoxon秩和检验
library(future)        # plan() 由 future 提供, future.apply 不重导出它(已实测)
library(future.apply)
# plan() 由 future 提供, 上一行的 future.apply 并不重导出它(已实测), 故已补 library(future)
plan(multisession)
pvalues <- future_lapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),group1)
  p=wilcox.test(gene~group1, data)$p.value
  return(p)
})
# 这里做了 FDR 校正 —— 优于本仓库多数脚本
fdr=p.adjust(pvalues,method = "fdr")

# 计算每个基因的倍数变化
group2=factor(group1)
conditionsLevel<-levels(group2)
# levels(factor(group1)) 按字典序排 -> [1]="0"(肿瘤), [2]="1"(正常), 与下面注释一致
dataCon1=count_norm[,c(which(group1==conditionsLevel[1]))] #肿瘤
dataCon2=count_norm[,c(which(group1==conditionsLevel[2]))] #正常
#肿瘤比正常，logFC大于0则为“基因在肿瘤上调”；若为正常比肿瘤，则logFC大于0表示“基因在正常中上调”
# logFC = log2(肿瘤均值 / 正常均值), 故 logFC>0 表示在【肿瘤】中更高,
# 与 14.DESeq_difference 及已统一方向后的 59.GEO_TCGA_Common_DEGs 一致
foldChanges=log2(rowMeans(dataCon1)/rowMeans(dataCon2)) 

# 基于FDR阈值的输出结果
pvalues0=t(as.data.frame(pvalues))
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues0, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
# ⚠ fdrThres 定义后未被使用, 下一行的阈值是写死的 0.05
fdrThres=0.05
#导出文件
write.table(outRst[with(outRst, ((log2foldChange> 1 | log2foldChange< (-1)) & FDR < 0.05 )),], file="wilcoxout.tsv",sep="\t", quote=F,row.names = T,col.names = T)
