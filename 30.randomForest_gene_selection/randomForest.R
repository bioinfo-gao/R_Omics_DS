# 随机森林特征基因筛选
#
# 输入: diffGeneExp.txt 行=基因, 列=样本(样本名需含 5 段以 - 分隔的字段, 见下)
# 输出: Random_Forest.pdf(误差曲线) / GeneIm.pdf(重要性) /
#       随机森林Genes.txt(入选基因) / imGeneExp.txt(入选基因表达)

# install.packages("randomForest")
#引用包
library(randomForest)

set.seed(123456)
setwd("C:/Users/zhen-/Code/R_code/R_Omics_DS/30.randomForest_gene_selection")

inputFile="diffGeneExp.txt"       #输入文件


#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
# ⚠ 分组取的是样本名第 5 段(\\5), 与本仓库其它脚本取第 4 段的约定不同。
#   若样本名不足 5 段, gsub 不匹配会原样返回整个样本名 ->
#   每个样本各成一类, randomForest 会把它当成 n 分类问题。不报错, 结果全无意义。
#   建议先 table(group) 确认只有两类再往下跑。
group=gsub("(.*)\\-(.*)\\-(.*)\\-(.*)\\-(.*)", "\\5", row.names(data))

#随机森林树

rf=randomForest(as.factor(group)~., data=data, ntree=1000)
# ⚠ 这里运行时写出的仍是中文名 森林.pdf; 仓库里同名文件已改名为 Random_Forest.pdf,
#   两者现在对不上。要保持一致需把本行也改成 Random_Forest.pdf。
pdf(file="森林.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()

#找出误差最小的点
# ⚠ OOB 误差曲线本身有随机波动, 取全局最小点作为 ntree 是在噪声上做选择;
#   通常只要 ntree 足够大到曲线走平即可, 不必挑最低点。
optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)

#查看基因的重要性
importance=importance(x=rf2)

#绘制基因的重要性图
pdf(file="GeneIm.pdf", width=6.2, height=5.8)
varImpPlot(rf2, main="")
dev.off()

#挑选疾病特征基因
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
# ⚠ MeanDecreaseGini > 2 是绝对阈值, 其数值随基因数与样本数变化, 不具可移植性。
#   下一行注释掉的「取前 30 名」是按秩选取, 换数据集时更稳定。
# 注: importance(rf2) 在默认 importance=FALSE 下是单列矩阵, order 后会降维成命名向量,
#     所以 names(rfGenes[...]) 才能取到基因名; 若改成 importance=TRUE 这行会失效。
rfGenes=names(rfGenes[rfGenes>2])     #挑选重要性评分大于2的基因
#rfGenes=names(rfGenes[1:30])         #挑选重要性评分最高的30个基因
write.table(rfGenes, file="随机森林Genes.txt", sep="\t", quote=F, col.names=F, row.names=F)

#输出重要基因的表达量
sigExp=t(data[,rfGenes])
sigExpOut=rbind(ID=colnames(sigExp),sigExp)
write.table(sigExpOut, file="imGeneExp.txt", sep="\t", quote=F, col.names=F)


