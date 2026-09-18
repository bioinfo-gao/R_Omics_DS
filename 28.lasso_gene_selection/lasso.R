# LASSO 逻辑回归筛选特征基因 (二分类: 由样本名第 5 段决定)
#
# 输入: diffGeneExp.txt 行=基因, 列=样本
# 输出: lambda.pdf / cvfit.pdf / LASSO.gene.txt / LASSO.geneExp.txt
#       LASSO.geneExp.txt 是 31.Machine_Learning_Modeling 的输入

setwd("C:/Users/zhen-/Code/R_code/R_Omics_DS/28.机器学习-lasso筛选基因")
set.seed(123)
library(glmnet)                   #???ð?
inputFile="diffGeneExp.txt"       #?????ļ?


rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=t(rt)

rt[1:2, 1:10]
x=as.matrix(rt)
# ⚠ 分组取样本名第 5 段(\\5)。若样本名不足 5 段, gsub 不匹配会原样返回整个样本名,
#   于是每个样本自成一类 —— glmnet 会报错或得到无意义结果。先 table(y) 确认只有两类。
y=gsub("(.*)\\-(.*)\\-(.*)\\-(.*)\\-(.*)", "\\5", row.names(rt))
fit=glmnet(x, y, family = "binomial", alpha=1)

pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

# 顶部已 set.seed(123), 故 cv.glmnet 的折划分可复现 —— 这一点优于本仓库其它几个建模脚本
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)

pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()


coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
# binomial 模型带截距, coef 的第 1 行是 (Intercept), 故用 [-1] 去掉。
# 注意 Cox 模型没有截距, 47/71 等脚本不需要这一步 —— 两者不能互抄。
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)

rt1=t(rt)
lassoexp=rt1[lassoGene,,drop=F]
lassoexp=as.data.frame(lassoexp)
colnames(lassoexp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(lassoexp))
write.table(lassoexp, file="LASSO.geneExp.txt", sep="\t", quote=F, row.names=T, col.names=T)

