# 单基因高低表达分组的无进展生存(PFS)曲线
#
# 输入: geneExp.txt(行=样本, 列含目标基因),
#       Survival_SupplementalTable_S1_xena_sp(Xena 的 TCGA 生存表, 取 PFI.time/PFI)
# 输出: <gene>.PFS.pdf

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("survival")
#install.packages("survminer")


library(limma)
library(survival)
library(survminer)

expFile="geneExp.txt"     #
cliFile="Survival_SupplementalTable_S1_xena_sp"      #?ٴ??????ļ?

getwd()
setwd("C:\\Users\\zhen-\\Code\\R_code\\R_For_DS_Omics\\05.PFS")      #???ù???Ŀ¼

#
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
# ⚠ 硬编码取第 3 列作为目标基因; 输入文件列序一变就换了基因且无提示
gene=colnames(rt)[3]
data=rt
# FIXME ⚠ 原写作 e=ifelse(...), 但下一行 cbind 用的是 Type ——
#   Type 在本脚本中从未定义, 运行必然报 object 'Type' not found。已改为 Type=。
Type=ifelse(data[,gene]>median(data[,gene]), "High", "Low")
data=cbind(as.data.frame(data), Type)

#??ȡ?ٴ??????ļ?
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[,c("PFI.time", "PFI")]
cli=na.omit(cli)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365
cli=as.matrix(cli)
# 把 Xena 的 barcode 截到患者层级, 以便与表达矩阵对齐
row.names(cli)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(cli))

#???ݺϲ???????????
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(as.data.frame(cli), data)

#?Ƚϸߵͱ?????֮???????????죬?õ???????pֵ
# 按 Type(高/低表达)比较 PFS; pValue 由卡方近似算出, 与 ggsurvplot 显示的是同一个
diff=survdiff(Surv(futime, fustat) ~ Type, data=rt)
pValue=1-pchisq(diff$chisq, df=1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=", sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Type, data = rt)

		
#????????????
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           surv.median.line = "hv",
		           legend.title=gene,
		           legend.labs=c("High level", "Low level"),
		           xlab="Time(years)",
		           ylab="Progression free survival",
		           break.time.by = 1,
		           palette=c("red", "blue"),
		           risk.table=F,
		       	   risk.table.title="",
		           risk.table.col = "strata",
		           risk.table.height=.25)

#????????????
pdf(file=paste0(gene, ".PFS.pdf"), width=5.5, height=5, onefile=FALSE)
print(surPlot)
dev.off()


