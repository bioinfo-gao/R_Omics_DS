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
# FIXME ⚠ 两个问题:
#   1) y=as.numeric(as.factor(group)) 把类别变成了 1/2 的【数值】, rfe 因此按【回归】处理,
#      这也是下面第 26 行取 RMSE 而不是 Accuracy 的原因。二分类应传 as.factor(group)。
#   2) methods="svmRadial" 多了一个 s。caretFuncs$fit 是 train(x, y, ...),
#      而 caret:::train.default 的 method 默认值是 "rf" —— 实测确认。
#      于是 methods= 被 ... 吞掉, 这里跑的其实是【随机森林】RFE, 不是 SVM。
#      正确写法是 method="svmRadial"。本文件未作改动, 由你确认后再改。
Profile=rfe(x=data,
            y=as.numeric(as.factor(group)),
            sizes = c(2,4,6,8, seq(10,40,by=3)),
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
            methods="svmRadial")

#绘制图形
pdf(file="SVM-RFE.pdf", width=6, height=5.5)
par(las=1)
x = Profile$results$Variables
y = Profile$results$RMSE
plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")
lines(x, y, col="darkgreen")
#标注交叉验证误差最小的点
wmin=which.min(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
points(wmin.x, wmin.y, col="blue", pch=16)
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)
dev.off()

#输出选择的基因
featureGenes=Profile$optVariables
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)
rt1=t(data)
# FIXME ⚠ lassoGene 在本脚本中从未定义 —— 它是 28.lasso_gene_selection/lasso.R 的变量。
#   全新会话里运行会报 object 'lassoGene' not found;
#   若先跑过 lasso.R 再跑本脚本, 则会静默拿 LASSO 选出的基因当成 SVM 的结果。
#   这里应当用第 38 行的 featureGenes。
svmexp=rt1[lassoGene,,drop=F]
svmexp=as.data.frame(svmexp)
colnames(svmexp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(svmexp))
# FIXME ⚠ 同上: 写出的是 lassoexp(LASSO 的表达矩阵), 而不是本脚本算出的 svmexp。
#   结果是 SVM.geneExp.txt 里装的其实是 LASSO 的结果。应改为 write.table(svmexp, ...)。
write.table(lassoexp, file="SVM.geneExp.txt", sep="\t", quote=F, row.names=T, col.names=T)
