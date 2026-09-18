#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("survival")
#install.packages("survminer")


setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/22.Immu_Survival")
library(limma)
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
Èread.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$time=cli$time/365

#???eSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(cli, data)
Type=Type[sameSample,,drop=F]

#???Tab=cbind(ID=row.names(rt), rt)
outTab=outTab[,-ncol(outTab)]
rt=outTab
sigGenes=c("time","event")
time=rt$time
data=rt[1]
exp=rt[,3:ncol(rt)]
rt=cbind(time,exp)
#???3 <- data.frame()
genes <- colnames(rt)[-c(1:2)]
plan(multisession)
system.time(res3 <- future_lapply(1:length(genes), function(i){
  group = ifelse(rt[,genes[i]]>median(rt[,genes[i]]),'high','low')
  if(length(table(group))==1) return(NULL)
  surv =as.formula(paste('Surv(time, event)~', 'group'))
  data = cbind(rt[,1:2],group)
  x = survdiff(surv, data = data)
  pValue=1-pchisq(x$chisq,df=1) 
  return(c(genes[i],pValue))
}))
res3 <- data.frame(do.call(rbind,res3))
names(res3 ) <- c('ID','pValue_log')
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


