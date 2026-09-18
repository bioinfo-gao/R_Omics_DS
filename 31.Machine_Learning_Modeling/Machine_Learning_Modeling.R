# 把特征基因表达二值化成「基因评分」, 再训练一个神经网络做二分类
#
# 输入: LASSO.geneExp.txt(来自 28.lasso_gene_selection), diff.txt(含 logFC 列)
# 输出: gene_score.txt, neural_network_model.pdf, neural.predict.txt
#
# ⚠⚠ 第 43 行有一处会让整个模型失去意义的跨脚本契约不匹配, 见该处说明。


library(limma)  
library(neuralnet)
library(NeuralNetTools)

expFile="LASSO.geneExp.txt"     #特征基因表达数据文件
diffFile="diff.txt"         #差异基因的文件


#读取表达文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#读取差异基因的文件
diffRT=read.table(diffFile, header=T, sep="\t", check.names=F, row.names=1)
# ⚠ 按行名对齐差异表; 若 diff.txt 缺少 data 里的某些基因, 对应行会变成 NA,
#   下面 logFC>0 / <0 的判断随之为 NA, data[NA,] 会产出整行 NA 且不报错。
#   稳妥做法: 先 intersect 行名, 或对齐后 stopifnot(!anyNA(diffRT$logFC))。
diffRT=diffRT[row.names(data),]

#基因评分
# 上调基因: 高于自身中位数记 1; 下调基因: 高于自身中位数记 0 —— 即按方向统一编码,
# 使「评分高」始终对应「更像处理组」
dataUp=data[diffRT[,"logFC"]>0,]
dataDown=data[diffRT[,"logFC"]<0,]
dataUp2=t(apply(dataUp,1,function(x)ifelse(x>median(x),1,0)))
dataDown2=t(apply(dataDown,1,function(x)ifelse(x>median(x),0,1)))

#输出基因评分的结果
outTab=rbind(dataUp2, dataDown2)
outTab=rbind(id=colnames(outTab), outTab)
write.table(outTab, file="gene_score.txt", sep="\t", quote=F, col.names=F)



#####构建模型
inputFile="gene_score.txt"       #输入文件
#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.data.frame(t(data))

#获取样品分组信息
# ⚠⚠ 跨脚本契约不匹配(已实测):
#   上游 28.lasso_gene_selection 在写 LASSO.geneExp.txt 时把样本名裁成了【4 段】
#   (gsub 保留 \\1-\\2-\\3-\\4), 例如 TCGA-AB-1234-01A;
#   而这里按【5 段】取第 5 段作为组别。段数不足 -> gsub 不匹配 -> 原样返回整个样本名。
#   后果: group 永远不等于 "con" 或 "treat", 下面两行使 con 与 treat 【全为 0】,
#   神经网络在全零标签上训练, 全程不报错。
#   未自动修改: 究竟哪一段承载组别标签, 取决于你的样本命名约定。
group=gsub("(.*)\\-(.*)\\-(.*)\\-(.*)\\-(.*)", "\\5", row.names(data))
data$con=ifelse(group=="con", 1, 0)
data$treat=ifelse(group=="treat", 1, 0)

#神经网络模型
# con+treat~. : 两个输出节点的多分类网络; "." 会自动排除已出现在左侧的 con/treat
fit=neuralnet(con+treat~., data, hidden=5)
fit$result.matrix
fit$weight
#plot(fit)

pdf(file="neural_network_model.pdf", width=10, height=10)
plotnet(fit)
dev.off()

#利用模型预测结果
# ⚠ compute() 传入了含 con/treat 两列的完整 data; 新版 neuralnet 推荐用 predict()。
#   另外若上面的分组失效, 下面 table(group, ...) 不是 2×2, 第 62-63 行的下标会越界报错。
net.predict=compute(fit, data)$net.result
net.prediction=c("con", "treat")[apply(net.predict, 1, which.max)]
predict.table=table(group, net.prediction)
predict.table
# 分别是对照组与处理组的预测准确率 —— 注意这是【训练集内】的准确率, 不是验证性能
conAccuracy=predict.table[1,1]/(predict.table[1,1]+predict.table[1,2])
treatAccuracy=predict.table[2,2]/(predict.table[2,1]+predict.table[2,2])
paste0("Con accuracy: ", sprintf("%.3f", conAccuracy))
paste0("Treat accuracy: ", sprintf("%.3f", treatAccuracy))

#输出预测结果
colnames(net.predict)=c("con", "treat")
outTab=rbind(id=colnames(net.predict), net.predict)
write.table(outTab, file="neural.predict.txt", sep="\t", quote=F, col.names=F)
