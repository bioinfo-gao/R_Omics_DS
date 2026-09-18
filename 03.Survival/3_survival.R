# 单基因生存筛选: KM + 单因素 Cox + 5 年生存率差值 三重过滤
#
# 输入: diffGeneExp.txt(行=基因, 列=样本), time.txt(行=样本, 含 futime/fustat)
# 输出: survival.xls(通过筛选的基因及 KM/HR/p), surSigExp.txt(这些基因的表达, 供 04.ROC 使用)
#
# ⚠ 本文件有一处在 R >= 4.2 上必然报错, 见第 67 行。

#install.packages("survival")
# if (!require("BiocManager", quietly = TRUE))
#         install.packages("BiocManager")
# BiocManager::install("limma")
# install.packages("tidyverse") # very slow, need 30 sec to finish 
# conda install conda-forge::r-tidyverse

library(tidyverse)
library(survival)
library(limma)
#setwd() 

expFile="diffGeneExp.txt"

cliFile="time.txt"  

rt2=read.table(expFile, header=T, check.names=F, row.names=1)

# 只取 -01A(原发肿瘤)列; 本脚本做的是肿瘤内部的预后分层, 不涉及癌旁
exp_data_T = rt2%>% dplyr::select(str_which(colnames(.), "-01A")) #

nT = ncol(exp_data_T) 

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$futime=cli$futime/365

# 下面 group1/conNum/treatNum/Type 算完之后【从未被使用】—— 是从差异分析脚本抄来的残留
group1=sapply(strsplit(colnames(rt2),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)
conNum=length(group1[group1==1])       
treatNum=length(group1[group1==0])     
Type=c(rep(1,conNum), rep(2,treatNum))


#ɾ��������Ʒ
tumorData=exp_data_T
tumorData=as.matrix(tumorData)
tumorData=t(tumorData)
# 把 TCGA barcode 截到前 3 段(患者层级), 以便与临床表按患者对齐
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)
#���ݺϲ���������
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt2=cbind(cli, data)
#����ϲ��������
outTab=cbind(ID=row.names(rt2), rt2)
# ⚠ 这一行丢掉的是【最后一列】。若 cli 恰好只有 futime/fustat 两列,
#   那么被丢掉的是最后一个基因, 且不会有任何提示。
#   若原意是去掉某个多余的临床列, 应按列名删除而不是按位置。
outTab=outTab[,-ncol(outTab)]


#�������
pFilter=0.05                                                                 #�����Թ�������
rt=outTab    #��ȡ�����ļ�                                                   #�������Ϊ��λ������30������Ϊ��λ������365
outTab=data.frame()
sigGenes=c("futime","fustat")
# 从第 4 列起是基因(1=ID, 2=futime, 3=fustat)
for(gene in colnames(rt[,4:ncol(rt)])){
# 表达几乎不变的基因直接跳过, 避免中位切分后两组完全相同
	   if(sd(rt[,gene])<0.1){
	     next}
# 以中位数二分高低表达; a 是逻辑向量, 下面 survdiff/survfit 直接按 ~a 分组
	   a=rt[,gene]<=median(rt[,gene])
	  
		 rt1=rt[a,]
	   b=setdiff(rownames(rt),rownames(rt1))
	   rt2=rt[b,]
	   surTab1=summary(survfit(Surv(futime, fustat) ~ 1, data = rt1))
	   surTab2=summary(survfit(Surv(futime, fustat) ~ 1, data = rt2))
	   survivalTab1=cbind(time=surTab1$time, surv=surTab1$surv,lower=surTab1$lower,upper=surTab1$upper)
	   survivalTab1=survivalTab1[survivalTab1[,"time"]<5,]
# FIXME ⚠ R >= 4.0 起 class(matrix) 返回 c("matrix","array") 两个元素,
#   于是 class(x)=="matrix" 是长度 2 的逻辑向量。R >= 4.2 会直接报
#   "the condition has length > 1"(在 R 4.5.2 上已实测), 循环无法跑完。
#   已改为 is.matrix()。
	   if(is.matrix(survivalTab1)){
	     survivalTab1=survivalTab1[nrow(survivalTab1),]
	   }
	   survivalTab2=cbind(time=surTab2$time, surv=surTab2$surv,lower=surTab2$lower,upper=surTab2$upper)
	   survivalTab2=survivalTab2[survivalTab2[,"time"]<5,]
	   if(is.matrix(survivalTab2)){
	     survivalTab2=survivalTab2[nrow(survivalTab2),]
	   }
# 取两组在 5 年内最后一个观测点的生存率之差; 若某组在 5 年内没有事件点会取到空值
	   fiveYearsDiff=abs(survivalTab1["surv"]-survivalTab2["surv"])

     #km����
# ⚠ ~a 里的 a 不是 rt 的列, 而是上面循环里的临时变量, 靠词法作用域找到。
#   写法可行但脆弱: 一旦把这段挪进函数, a 就可能解析到别的东西。
	   diff=survdiff(Surv(futime, fustat) ~a,data = rt)
	   pValue=1-pchisq(diff$chisq,df=1)
	   fit=survfit(Surv(futime, fustat) ~ a, data = rt)
	   #cox����
	   cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
	   coxSummary = summary(cox)
	   coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	
# 三重过滤: KM 的 p、Cox 的 p、以及 5 年生存率差值 > 0.15。
# ⚠ 三个条件都用未校正的原始 p, 基因数一多假阳性会很可观。
	   if((pValue<pFilter) & (coxP<pFilter) & (fiveYearsDiff>0.15)){
	       sigGenes=c(sigGenes,gene)
	       outTab=rbind(outTab, 
	                    cbind(gene=gene,
	                          KM=pValue,
	                          HR=coxSummary$conf.int[,"exp(coef)"],
	                          HR.95L=coxSummary$conf.int[,"lower .95"],
	                          HR.95H=coxSummary$conf.int[,"upper .95"],
			                      coxPvalue=coxP) )
		 }
}
write.table(outTab,file="survival.xls",sep="\t",row.names=F,quote=F)    #��������pֵ�����ļ�
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="surSigExp.txt",sep="\t",row.names=F,quote=F)

