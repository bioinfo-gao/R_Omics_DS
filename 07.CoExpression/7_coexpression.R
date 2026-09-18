# 目标基因的全基因组共表达分析 (Pearson)
#
# 输入: combined_RNAseq_FPKM.txt; 顶部 gene 变量指定目标基因
# 输出: corResult.txt(全部结果), corSig.txt(达到阈值者), 以及每个达标基因一张 cor.<gene>.pdf
#
# ⚠⚠ 本文件原状态【完全无法运行】, 有三处阻断性问题, 均已修复, 见下方说明。

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")


#引用包
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)

gene="ARL15"               #目标基因的名称
corFilter=0.3            #相关系数的过滤条件
pFilter=0.05             #相关性检验pvalue的过滤条件
expFile="combined_RNAseq_FPKM.txt"      #表达数据文件


#读取输入文件，并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
# FPKM 上的弱过滤; 注意下面还会再取 log2
data=data[rowMeans(data)>1,]

#删掉正常样品
# 只保留肿瘤样本(barcode 第 4 段首字符为 0)
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
# log2(x+1) 之后再算 Pearson —— 相关性因此是在对数尺度上定义的
data=log2(data+1)

#提取目标基因表达量
x=as.numeric(data[gene,])
#对基因进行循环，进行相关性检验
# FIXME ⚠ 下面这段并行化尝试是坏的, 已整段注释掉。它有四处问题:
#   1) rt2 在本脚本中从未定义;
#   2) plan()/future_lapply() 来自 future 包, 而本文件没有 library(future);
#   3) 括号位置错了 —— future_lapply(rownames(data)) 只传了一个参数,
#      匿名函数实际是作为第二个参数传给了 system.time;
#   4) next 只能用在循环里, 放在函数体内会报错。
#   真正在工作的是下面那个普通 for 循环。
# res3 <- data.frame()
# genes <- colnames(rt2)[-c(1:2)]
# plan(multisession)
# system.time(res3 <- future_lapply(rownames(data)), function(j){
#   if(gene==j){next}
#   y=as.numeric(data[j,])
#   corT=cor.test(x, y, method = 'pearson')
#   cor=corT$estimate
#   pvalue=corT$p.value
#   outTab=rbind(outTab, cbind(Query=gene, Gene=j, cor, pvalue))
#   #保存满足条件的基因
#   if((abs(cor)>corFilter) & (pvalue<pFilter)){
#     #可视化
#     df1=as.data.frame(cbind(x,y))
#     p1=ggplot(df1, aes(x, y)) + 
#       xlab(paste0(gene, " expression"))+ ylab(paste0(j, " expression"))+
#       geom_point()+ geom_smooth(method="lm", formula=y~x) + theme_bw()+
#       stat_cor(method = 'pearson', aes(x =x, y =y))
#     pdf(file=paste0("cor.", j, ".pdf"), width=5, height=4.6)
#     print(p1)
#     dev.off()
#   }
# })
# res3 <- data.frame(do.call(rbind,res3))
# 逐基因做相关性检验; 达到阈值的顺便出一张散点图
# ⚠ 未做多重检验校正: 全转录组两两相关时 p<0.05 会产生大量假阳性,
#   若要下结论请改用 BH 校正后的 q 值。
outTab=data.frame()   # 必须先初始化, 否则 rbind 报 object 'outTab' not found
for(j in rownames(data)){
	if(gene==j){next}
    y=as.numeric(data[j,])
	corT=cor.test(x, y, method = 'pearson')
	cor=corT$estimate
	pvalue=corT$p.value
	outTab=rbind(outTab, cbind(Query=gene, Gene=j, cor, pvalue))
	#保存满足条件的基因
	if((abs(cor)>corFilter) & (pvalue<pFilter)){
		#可视化
		df1=as.data.frame(cbind(x,y))
		p1=ggplot(df1, aes(x, y)) + 
			xlab(paste0(gene, " expression"))+ ylab(paste0(j, " expression"))+
			geom_point()+ geom_smooth(method="lm", formula=y~x) + theme_bw()+
			stat_cor(method = 'pearson', aes(x =x, y =y))
		pdf(file=paste0("cor.", j, ".pdf"), width=5, height=4.6)
		print(p1)
		dev.off()
	}
}

#输出相关性结果文件
write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)
# ⚠ 这一行需要 outTab 是 data.frame 才能用 $ 取列。
#   原代码没有初始化 outTab, rbind 会把结果拼成【字符矩阵】, $cor 会直接报
#   "$ operator is invalid for atomic vectors"(已实测)。
#   现已在循环前加 outTab=data.frame() 解决。
outTab=outTab[abs(as.numeric(outTab$cor))>corFilter & as.numeric(outTab$pvalue)<pFilter,]
write.table(file="corSig.txt", outTab, sep="\t", quote=F, row.names=F)



