# 合并 TCGA 与 GTEx 表达矩阵并用 ComBat 去批次
#
# 输入: combined_RNAseq_TPM_log.txt(TCGA), GTExNormalExp.txt(GTEx 正常组织)
# 输出: TCGA_GTEx_normalize.txt, boxbind.pdf, boxbind_normal.pdf
#
# ⚠⚠ 第 35 行涉及一个 TCGA+GTEx 合并的经典陷阱, 务必先读该处说明再用本脚本的结果。



#引用包
library(limma)
library(sva)

#读取基因表达文件,并对数据进行处理
# ⚠ 注意这一侧【没有】像下面 GTEx 那样显式设置行名。
#   它能工作的前提是: 该文件由 write.table(row.names=TRUE, col.names=TRUE) 写出,
#   表头字段比数据行少一个, read.table 因此自动把第一列当作行名。
#   若文件不是这种「表头错位」格式, rownames 会变成 1,2,3..., 第 19 行的 intersect
#   将得到空集, TCGA 侧变成 0 行 —— 不报错, 但后面全是空的。
#   稳妥做法: 改成与下面 GTEx 一致的显式写法(rownames(rt)=rt[,1] 那四行)。
rt=read.table("combined_RNAseq_TPM_log.txt", header=T, sep="\t", check.names=F)
data=as.data.frame(rt)
data0=avereps(data)
GTEx=read.table("GTExNormalExp.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(GTEx)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data2=as.data.frame(data)
# 只保留两套数据共有的基因
sameGene=intersect(rownames(data0),rownames(data2))
TCGA=data0[sameGene,,drop=F]
nTCGA = ncol(TCGA) 
GTEx=data2[sameGene,,drop=F]
nGTEx = ncol(GTEx) 
data0=cbind(TCGA,GTEx)
pdf(file = "boxbind.pdf")
# ⚠ 输入文件名已带 _log, 若其确为 log 尺度, 这里的 log2(x+1) 就是第二次取对数
boxplot(log2(data0+1),outline=T, notch=T,las=2)
dev.off()
TCGA_smaple=as.data.frame(colnames(TCGA))
colnames(TCGA_smaple)="group"
TCGA_smaple$batch="1"
GTEx_sample=as.data.frame(colnames(GTEx))
colnames(GTEx_sample)="group"
GTEx_sample$batch="2"
sample=rbind(TCGA_smaple,GTEx_sample)
# ⚠⚠ 关键: mod 只含截距, 意味着告诉 ComBat「除批次外没有需要保留的生物学变量」。
#   但本脚本合并的是【TCGA 肿瘤】与【GTEx 正常组织】——
#   批次(TCGA/GTEx)与生物学分组(肿瘤/正常)是【完全共线】的。
#   在这种情形下去批次会把肿瘤与正常的真实差异一并抹掉, 之后再做差异分析
#   往往得到偏小甚至方向错误的结果, 而过程不会有任何报错。
#   若确实要合并这两个来源, 应把生物学分组放进 mod(如 model.matrix(~group)),
#   并清楚这在共线情况下仍然是有风险的。
#   未自动修改: 这属于研究设计层面的决定。
mod = model.matrix(~1,data = sample)
batch=sample$batch
# par.prior=T 用参数化经验贝叶斯估计批次参数(样本量较大时更稳且更快)
data0_combat=ComBat(dat=data0,batch = batch,mod=mod,par.prior = T)
write.table(data0_combat,file = "TCGA_GTEx_normalize.txt",sep="\t",quote=F)
pdf(file = "boxbind_normal.pdf")
boxplot(log2(data0_combat+1),outline=T, notch=T,las=2)
dev.off()
