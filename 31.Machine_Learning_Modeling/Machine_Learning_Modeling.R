
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
diffRT=diffRT[row.names(data),]

#基因评分
dataUp=data[diffRT[,"logFC"]>0,]
dataDown=data[diffRT[,"logFC"]<0,]
dataUp2=t(apply(dataUp,1,function(x)ifelse(x>median(x),1,0)))
dataDown2=t(apply(dataDown,1,function(x)ifelse(x>median(x),0,1)))

#输出基因评分的结果
outTab=rbind(dataUp2, dataDown2)
outTab=rbind(id=colnames(outTab), outTab)
write.table(outTab, file="基因评分.txt", sep="\t", quote=F, col.names=F)



#####构建模型
inputFile="基因评分.txt"       #输入文件
#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.data.frame(t(data))

#获取样品分组信息
group=gsub("(.*)\\-(.*)\\-(.*)\\-(.*)\\-(.*)", "\\5", row.names(data))
data$con=ifelse(group=="con", 1, 0)
data$treat=ifelse(group=="treat", 1, 0)

#神经网络模型
fit=neuralnet(con+treat~., data, hidden=5)
fit$result.matrix
fit$weight
#plot(fit)

pdf(file="神经网络模型.pdf", width=10, height=10)
plotnet(fit)
dev.off()

#利用模型预测结果
net.predict=compute(fit, data)$net.result
net.prediction=c("con", "treat")[apply(net.predict, 1, which.max)]
predict.table=table(group, net.prediction)
predict.table
conAccuracy=predict.table[1,1]/(predict.table[1,1]+predict.table[1,2])
treatAccuracy=predict.table[2,2]/(predict.table[2,1]+predict.table[2,2])
paste0("Con accuracy: ", sprintf("%.3f", conAccuracy))
paste0("Treat accuracy: ", sprintf("%.3f", treatAccuracy))

#输出预测结果
colnames(net.predict)=c("con", "treat")
outTab=rbind(id=colnames(net.predict), net.predict)
write.table(outTab, file="neural.predict.txt", sep="\t", quote=F, col.names=F)
