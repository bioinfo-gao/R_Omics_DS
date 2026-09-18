# 多指标 ROC 对比: riskScore 与各临床因素在同一张图上比 1 年预测能力
#
# 输入: Risk.txt(同 54), Clinical.txt
# 输出: ROC.pdf



library(survivalROC)

# ⚠ 形参默认值写的是 null(小写), R 无此对象; 因调用时都显式传值而未触发。已改为 NULL。
bioROC=function(riskFile=NULL,cliFile=NULL,outFile=NULL){
		risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        #读取风险文件
		cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          #读取临床文件
		sameSample=intersect(row.names(cli),row.names(risk))
		risk=risk[sameSample,]
		cli=cli[sameSample,]
# ⚠ 与 54 相同的按位置取列, 依赖 Risk.txt 的列序
		rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
		rocCol=rainbow(ncol(rt)-2)
		aucText=c()
		
		#绘制risk score的ROC曲线
		pdf(file=outFile,width=6,height=6)
		par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
# predict.time=1 表示 1 年; 这要求 futime 的单位是【年】(上游 46/71 已除以 365)
		roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")
		plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
		  xlab="False positive rate", ylab="True positive rate",
		  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
		aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")"))
		abline(0,1)
		
		#绘制其他临床性状的ROC曲线
		j=1
# 遍历临床变量(排除末列 riskScore, 它已单独画过)。
# ⚠ survivalROC 的 marker 必须是数值: 临床变量若是 "Stage I" 这类字符会直接报错,
#   需要先自行编码成有序数值。
		for(i in colnames(rt[,3:(ncol(rt)-1)])){
			roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =1, method="KM")
			j=j+1
			aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
			lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
		}
# 图例按 rocCol 顺序排列, 与曲线绘制顺序一致
		legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
		dev.off()
}
bioROC(riskFile="Risk.txt",cliFile="Clinical.txt",outFile="ROC.pdf")

