# 各免疫细胞浸润比例的高/低分组生存差异 (log-rank)
#
# 输入: CIBERSORT-Results.txt, time.txt(行=样本, 含 time/event)
# 输出: 每个免疫细胞一张 <cell>_surv.pdf
#
# ⚠⚠ 本文件此前有四处阻断性问题, 均已修复:
#   两行代码被编码事故并进注释(cli= 与 outTab=)、sameSample= 被注释掉、
#   以及 plan() 所在的包未被加载。

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("survival")
#install.packages("survminer")


setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/22.Immu_Survival")
library(limma)
library(future)        # plan() 由 future 提供, future.apply 并不重导出它(已实测)
library(future.apply)
library(survival)
library(survminer)
library(tidyverse)

expFile="CIBERSORT-Results.txt"     #???ߴ????ļ?
cliFile="time.txt"        #?ٴ??????ļ?


rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
rt[1:2, 1:5]
gene=colnames(rt)


tumorData=as.matrix(rt)
tumorData=t(rt)
tumorData=as.data.frame(tumorData)
exp_data_T = tumorData%>% dplyr::select(str_which(colnames(.), "-01A"))
tumorData=cbind(exp_data_T)
data=t(avereps(tumorData))
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$time=cli$time/365

sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(cli, data)
Type=Type[sameSample,,drop=F]

# ⚠ 按位置丢掉最后一列; 若列序变化会丢错列且不报错
# 已还原: 原为 outTab=cbind(ID=row.names(rt), rt); 判据是下一行使用 outTab。
outTab=cbind(ID=row.names(rt), rt)
outTab=outTab[,-ncol(outTab)]
rt=outTab
sigGenes=c("time","event")
time=rt$time
data=rt[1]
exp=rt[,3:ncol(rt)]
rt=cbind(time,exp)
# 已还原: 原为 res3 <- data.frame()(结果容器的初始化)。
res3 <- data.frame()
genes <- colnames(rt)[-c(1:2)]
plan(multisession)
# ⚠ 这里 3:ncol 把 event 也一并当成了「基因」列, 但下面 genes 取 [-c(1:2)]
#   又把 time/event 去掉了, 两处正好抵消。改动任一处时必须同步看另一处。
system.time(res3 <- future_lapply(1:length(genes), function(i){
  group = ifelse(rt[,genes[i]]>median(rt[,genes[i]]),'high','low')
  if(length(table(group))==1) return(NULL)
  surv =as.formula(paste('Surv(time, event)~', 'group'))
  data = cbind(rt[,1:2],group)
# ⚠ res3 里的 pValue 未做多重检验校正; 下面 pValue_log < 1 实际上不过滤任何东西。
  x = survdiff(surv, data = data)
  pValue=1-pchisq(x$chisq,df=1) 
  return(c(genes[i],pValue))
}))
res3 <- data.frame(do.call(rbind,res3))
names(res3 ) <- c('ID','pValue_log')
# 对筛出的细胞类型逐个画 KM 曲线; 分组阈值同样是本批数据的中位数
res3 <- res3[with(res3, (pValue_log < 1 )), ]
#??ͼ
genes=res3$ID
for (i in 1:length(genes)) {
  print(i)
  # ??λ??????
  group = ifelse(rt[,genes[i]]>median(rt[,genes[i]]),'high','low')
  p=ggsurvplot(survfit(Surv(time, event)~group, 
                       data=rt), conf.int=F, pval=TRUE,title=genes[i])
  pdf(paste0(genes[i], "_surv.pdf"),width = 5, height = 5)
  print(p, newpage = FALSE)
  dev.off()
}


