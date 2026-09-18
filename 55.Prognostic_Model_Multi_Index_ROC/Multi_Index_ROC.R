

library(survivalROC)

bioROC=function(riskFile=null,cliFile=null,outFile=null){
		risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        #读取风险文件
		cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          #读取临床文件
		sameSample=intersect(row.names(cli),row.names(risk))
		risk=risk[sameSample,]
		cli=cli[sameSample,]
		rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
		rocCol=rainbow(ncol(rt)-2)
		aucText=c()
		
		#绘制risk score的ROC曲线
		pdf(file=outFile,width=6,height=6)
		par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
		roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")
		plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
		  xlab="False positive rate", ylab="True positive rate",
		  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
		aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")"))
		abline(0,1)
		
		#绘制其他临床性状的ROC曲线
		j=1
		for(i in colnames(rt[,3:(ncol(rt)-1)])){
			roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =1, method="KM")
			j=j+1
			aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
			lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
		}
		legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
		dev.off()
}
bioROC(riskFile="Risk.txt",cliFile="Clinical.txt",outFile="ROC.pdf")

