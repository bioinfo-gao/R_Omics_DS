# 单基因诊断 ROC 曲线 (肿瘤 vs 癌旁)
#
# 输入: combined_RNAseq_counts.txt (原始 counts) + 顶部 gene 变量指定的基因
# 输出: ROC.<gene>.pdf




library(pROC)
library(limma)
library(tidyverse)

#输入目标基因
gene="THBS2"
#读取输入文件
# ⚠ 这里读的是原始 counts, 未做 CPM/TPM 归一化, 样本间测序深度差异未校正。
#   单基因 ROC 对深度差异敏感, 建议改用 TPM 或 CPM 输入。
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
# 从 barcode 第 4 段首字符取样本类型: 0=肿瘤, 1=正常(2=复发并入 1)
group=sapply(strsplit(rownames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
group=as.data.frame(group)  
data$group=group$group
y=data$group


#绘制ROC曲线
# ⚠ 两个默认行为会静默影响结论:
#   1) y 是字符型 "0"/"1", pROC 按字母序取 "0"(肿瘤)为对照、"1"(正常)为病例,
#      方向与直觉相反 —— AUC 描述的是「用该基因把正常判为阳性」的能力。
#   2) roc() 默认 direction="auto", 会自动翻转比较方向以保证 AUC>=0.5,
#      因此 AUC 永远不会低于 0.5, 不能据此说明基因有判别力。
#      要如实反映方向请显式写 direction="<" 或 ">"。
roc1=roc(y, as.numeric(data[,gene]))
# bootstrap 置信区间已固定种子, 可复现
set.seed(123)   # bootstrap 抽样固定, 保证 95% CI 可复现
ci1=ci.auc(roc1, method="bootstrap")
ciVec=as.numeric(ci1)
pdf(file=paste0("ROC.",gene,".pdf"), width=5, height=5)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=gene)
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()

