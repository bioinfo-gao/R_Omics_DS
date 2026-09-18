# 机器学习模型预测结果的 ROC 曲线
#
# 输入: ML_prediction.txt —— 本仓库内没有任何脚本生成该文件, 需自行准备。
#       ⚠ 原文件名为中文「机器学习预测结果.txt」, 已改为 ASCII;
#         若你手上已有该数据文件, 需同步改名后本脚本才能读到。
# 输出: ROC.pdf



library(pROC)                   
inputFile="ML_prediction.txt"      



rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
# ⚠ 两层假设叠在一起, 任一不成立都会静默出错:
#   1) 样本名必须是 5 段短横线格式, 否则 gsub 原样返回整个名字;
#   2) 第 5 段必须恰好是 "con", 否则 ifelse 会把【所有】样本判为 1,
#      y 只剩一个水平, roc() 会报 'No control observation'。
y=gsub("(.*)\\-(.*)\\-(.*)\\-(.*)\\-(.*)", "\\5", row.names(rt))
y=ifelse(y=="con", 0, 1)

#绘制ROC曲线
# ⚠ 三点:
#   1) 硬编码取第 2 列作为预测值, 输入文件列序一变就取错列;
#   2) roc() 默认 direction="auto" 会自动翻转方向以保证 AUC>=0.5,
#      所以 AUC 永远不低于 0.5, 不能据此说明模型有判别力;
#   3) 图标题写的是 Train group —— 这是训练集内部的 ROC, 不是验证性能。
roc1=roc(y, as.numeric(rt[,2]))
# bootstrap 置信区间已固定种子, 可复现
set.seed(123)   # bootstrap 抽样固定, 保证 95% CI 可复现
ci1=ci.auc(roc1, method="bootstrap")
ciVec=as.numeric(ci1)
pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main="Train group")
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()


