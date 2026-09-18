#读取TPM矩阵
readCount<-read.table(file="combined_RNAseq_TPM.txt", header = T, row.names = 1, stringsAsFactors = F,check.names = F)
# edger 标准化并删除低表达基因
library(edgeR)
library(limma)
group1=sapply(strsplit(colnames(readCount),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)
y <- DGEList(counts=readCount,group=group1)
#删除过低表达量基因
keep <- filterByExpr(y)
y <- y[keep,keep.lib.sizes=FALSE]
##进行TMM标准化并转移到CPM（百万计数）
y <- calcNormFactors(y,method="TMM")
count_norm=cpm(y)
count_norm<-as.data.frame(count_norm)

# 对每个基因进行Wilcoxon秩和检验
library(future.apply)
plan(multisession)
pvalues <- future_lapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),group1)
  p=wilcox.test(gene~group1, data)$p.value
  return(p)
})
fdr=p.adjust(pvalues,method = "fdr")

# 计算每个基因的倍数变化
group2=factor(group1)
conditionsLevel<-levels(group2)
dataCon1=count_norm[,c(which(group1==conditionsLevel[1]))] #肿瘤
dataCon2=count_norm[,c(which(group1==conditionsLevel[2]))] #正常
#肿瘤比正常，logFC大于0则为“基因在肿瘤上调”；若为正常比肿瘤，则logFC大于0表示“基因在正常中上调”
foldChanges=log2(rowMeans(dataCon1)/rowMeans(dataCon2)) 

# 基于FDR阈值的输出结果
pvalues0=t(as.data.frame(pvalues))
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues0, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
fdrThres=0.05
#导出文件
write.table(outRst[with(outRst, ((log2foldChange> 1 | log2foldChange< (-1)) & FDR < 0.05 )),], file="wilcoxout.tsv",sep="\t", quote=F,row.names = T,col.names = T)
