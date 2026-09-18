# LASSO-Cox 预后模型: 由单因素筛选结果建模并输出风险评分
#
# 输入: UniSigExp.txt (来自 44.Prognosis_Related_Gene_Selection 的单因素 Cox 筛选)
#       列序契约: 1=样本id(→rownames) 2=futime 3=fustat 4+=基因表达
# 输出: genes_Coef.txt (入选基因及其系数), Risk.txt (id/futime/fustat/基因/riskScore/risk)
#       Risk.txt 是 54.Independent_Prognostic_Analysis 与 55.*_Multi_Index_ROC 的输入


library("glmnet")
library("survival")

set.seed(123)   # cv.glmnet 的折划分是随机的; 不设种子则每次入选基因都可能不同

###设置工作目录
# ⚠ setwd("") 会直接报错(cannot change working directory), 运行前必须填上真实路径
# setwd("")   # 原为空字符串会直接报错; 填入真实路径后再取消注释    

###读取文件
rt=read.table("UniSigExp.txt",header=T,sep="\t",row.names=1)            
rt$futime=rt$futime/365

###构建模型
# 3:ncol 按位置取基因列, 依赖 futime/fustat 恰在前两列
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 1000)
# 折划分已由文件顶部的 set.seed(123) 固定, 结果可复现。
# ⚠ maxit=1000 远低于 glmnet 默认的 1e5, 可能未收敛就停止(只给 warning 不报错)。
cvfit=cv.glmnet(x, y, family="cox", maxit = 1000)

###输出相关基因系数
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="genes_Coef.txt",sep="\t",quote=F,row.names=F)

###输出风险值
# lassoGene 与 actCoef 由同一 index 取出, 列序对齐; 换写法时务必保持这一点
trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
# ⚠ 中位切点来自本批样本, 且 riskScore 由建模同一批数据算出(in-sample),
#   high/low 的分离程度天然乐观, 不能作为模型验证结果。
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="Risk.txt",sep="\t",quote=F,row.names=F)

