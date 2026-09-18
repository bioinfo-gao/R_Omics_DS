# LASSO-Cox 建模并输出风险评分
#
# 输入: UniSigExp.txt (单因素 Cox 筛选结果; 列 1=样本id→rownames 2=futime 3=fustat 4+=基因)
# 输出: geneCoef.txt(入选基因与系数), Risk.txt(id/futime/fustat/基因/riskScore/risk)
#
# 注: 本文件与 71.Build_Lasso_Prognostic_Model 几乎相同, 差别是那边多了 setwd 与 maxit,
#     且输出文件名为 genes_Coef.txt。两处是同一套流程的两个副本, 改动时注意同步。



library("glmnet")
library("survival")

set.seed(123)   # cv.glmnet 的折划分是随机的; 不设种子则每次入选基因都可能不同

rt=read.table("UniSigExp.txt",header=T,sep="\t",row.names=1)            
rt$futime=rt$futime/365
# 3:ncol 按位置取基因列, 依赖 futime/fustat 恰在前两列

#构建模型
x=as.matrix(rt[,c(3:ncol(rt))])
# 折划分随机性已由顶部 set.seed(123) 固定; 不设种子时每次入选基因都可能不同
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox")
cvfit=cv.glmnet(x, y, family="cox")

#输出相关基因系数
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)
# lassoGene 与 actCoef 取自同一 index, 列序对齐

#输出风险值
trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
# ⚠ 中位切点与 riskScore 都来自建模同一批样本(in-sample), 分离程度天然乐观,
#   不能当作模型验证结果。
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="Risk.txt",sep="\t",quote=F,row.names=F)
