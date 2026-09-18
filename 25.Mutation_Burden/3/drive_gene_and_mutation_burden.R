# 驱动基因表达 与 肿瘤突变负荷(TMB) 的关系
#
# 输入: geneExp.txt(目标基因表达), TMB.txt(每个样本的突变负荷)
# 输出: 表达高低分组与 TMB 的比较图
#
# ⚠ 第 20/25 行原有尾部反斜杠 bug, 已修复(见下)。

setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/25.肿瘤突变负荷分析/3")


#???ð?
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)

expFile="geneExp.txt"     #?????????ļ?
tmbFile="TMB.txt"         #????ͻ???????ļ?

#??ȡ?????????ļ?
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[1]

#ɾ????????Ʒ
tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
# ⚠ 已修复: 两处 gsub 的替换串结尾都多一个反斜杠, 会给样本名追加字面反斜杠。
#   本文件里因为【两侧都被同样地加了反斜杠】, intersect 仍能匹配上, 侥幸未出错;
#   但只要有一侧改用正确写法, 匹配就会全部落空。现已一并改正。
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

#??ȡ????ͻ?为?ɵ??ļ?
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(tmb)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tmb))
tmb=avereps(tmb)
#???ݺϲ???????????
sameSample=intersect(row.names(data), row.names(tmb))
data=data[sameSample,,drop=F]
tmb=tmb[sameSample,,drop=F]
rt=cbind(data, tmb)

#?????Է???
x=as.numeric(rt[,gene])
y=log2(as.numeric(rt[,"total_perMB_log"])+1)
df1=as.data.frame(cbind(x,y))
corT=cor.test(x, y, method="spearman")
p1=ggplot(df1, aes(x, y)) + 
			xlab(paste0(gene, " expression"))+ylab("Tumor mutation burden")+
			geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
			stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))

#??????????ͼ??
pdf(file="cor.pdf",width=5,height=5)
print(p2)
dev.off()



