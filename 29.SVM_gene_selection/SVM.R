# SVM-RFE 特征基因筛选
#
# 输入: diffGeneExp.txt 行=基因, 列=样本
# 输出: SVM-RFE.pdf / SVM-RFE.gene.txt / SVM.geneExp.txt
#
# ⚠⚠ 本文件有三处问题, 其中两处会让输出内容根本不是 SVM 的结果。见下方 FIXME。


#引用包
library(e1071)
library(kernlab)
library(caret)

set.seed(123)
inputFile="diffGeneExp.txt"        #输入文件

#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
# ⚠ 同 28/30: 取样本名第 5 段, 不足 5 段时会静默退化成「每样本一类」
group=gsub("(.*)\\-(.*)\\-(.*)\\-(.*)\\-(.*)", "\\5", row.names(data))

#SVM-RFE分析
# 三处已修复(原状态会让本脚本既不是 SVM 也不是分类):
#   1) y 原为 as.numeric(as.factor(group)) 的 1/2 数值, rfe 按回归处理 -> 改为 factor;
#   2) methods="svmRadial" 多一个 s, 被 ... 吞掉, caret:::train.default 默认 method="rf"
#      (已实测确认) -> 改为 method=;
#   3) 因而下面的选优指标由 RMSE 最小改为 Accuracy 最大。
Profile=rfe(x=data,
            y=as.factor(group),
            sizes = c(2,4,6,8, seq(10,40,by=3)),
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
            method="svmRadial")

#绘制图形
pdf(file="SVM-RFE.pdf", width=6, height=5.5)
par(las=1)
x = Profile$results$Variables
y = Profile$results$Accuracy
plot(x, y, xlab="Variables", ylab="Accuracy (Cross-Validation)", col="darkgreen")
lines(x, y, col="darkgreen")
#标注交叉验证准确率最高的点(分类任务取最大值)
wmin=which.max(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
points(wmin.x, wmin.y, col="blue", pch=16)
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)
dev.off()

#输出选择的基因
featureGenes=Profile$optVariables
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)
rt1=t(data)
# 已修复: 原先写的是 lassoGene(28.lasso_gene_selection 的变量, 本脚本从未定义)。
# 同一会话里先跑过 lasso.R 的话, 会静默把 LASSO 的基因当成 SVM 的结果。
svmexp=rt1[featureGenes,,drop=F]
svmexp=as.data.frame(svmexp)
colnames(svmexp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(svmexp))
# 已修复: 原先写出的是 lassoexp(LASSO 的表达矩阵), 导致 SVM.geneExp.txt 装的是 LASSO 结果。
write.table(svmexp, file="SVM.geneExp.txt", sep="\t", quote=F, row.names=T, col.names=T)
