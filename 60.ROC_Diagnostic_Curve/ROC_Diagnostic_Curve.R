


library(pROC)
library(limma)
library(tidyverse)

#输入目标基因
gene="THBS2"
#读取输入文件
rt=read.table(file ="combined_RNAseq_counts.txt" , header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=as.data.frame(data)
#以01A和11A分组，正常放前面，肿瘤放后面
exp_data_T = data%>% dplyr::select(str_which(colnames(.), "-01A")) # 匹配列名或用下示写法
nT = ncol(exp_data_T) 
exp_data_N = data%>% dplyr::select(str_which(colnames(.), "-11A"))
nN = ncol(exp_data_N) 
data= cbind(exp_data_N, exp_data_T)
data=avereps(data)
data=t(data[gene,,drop=F])
data=as.data.frame(data)
group=sapply(strsplit(rownames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
group=as.data.frame(group)  
data$group=group$group
y=data$group


#绘制ROC曲线
roc1=roc(y, as.numeric(data[,gene]))
ci1=ci.auc(roc1, method="bootstrap")
ciVec=as.numeric(ci1)
pdf(file=paste0("ROC.",gene,".pdf"), width=5, height=5)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=gene)
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()

