# 独立预后分析: 判断 riskScore 在校正临床因素后是否仍然独立预测生存
#
# 输入: Risk.txt(来自 46 或 71; 列序 1=futime 2=fustat ... 倒数第2=riskScore 末列=risk),
#       Clinical.txt(行=样本, 列=临床变量)
# 输出: unCox.txt / muCox.txt 两张表, unForest.pdf / muForest.pdf 两张森林图


library(survival)
risk=read.table("Risk.txt",header=T,sep="\t",check.names=F,row.names=1)        #读取风险文件
cli=read.table("Clinical.txt",sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
# ⚠ 完全按位置取列: 1=futime, 2=fustat, 倒数第 2 = riskScore。
#   Risk.txt 的列序一旦变化, 这里会静默取到别的列。
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])

#单因素独立预后分析
uniTab=data.frame()
# 单因素: 每个临床变量和 riskScore 各自单独进模型
for(i in colnames(rt[,3:ncol(rt)])){
	 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	 coxSummary = summary(cox)
	 uniTab=rbind(uniTab,
	              cbind(id=i,
	              HR=coxSummary$conf.int[,"exp(coef)"],
	              HR.95L=coxSummary$conf.int[,"lower .95"],
	              HR.95H=coxSummary$conf.int[,"upper .95"],
	              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
	              )
}
write.table(uniTab,file="unCox.txt",sep="\t",row.names=F,quote=F)

#多因素独立预后分析
# 多因素: "." 把 rt 里除 futime/fustat 外的所有列一并放入, 不做任何变量选择。
# ⚠ 临床变量若是字符型会被当作因子展开成多个哑变量; 缺失值会导致该样本被整行丢弃。
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="muCox.txt",sep="\t",row.names=F,quote=F)


############绘制森林图函数############
# ⚠ 形参默认值写的是 null(小写), R 里没有这个对象, 正确写法是 NULL。
#   因为三个参数在调用时都显式传了值, 默认值从未被求值, 所以侥幸不报错。已改为 NULL。
bioForest=function(coxFile=NULL,forestFile=NULL,forestCol=NULL){
		rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
		gene <- rownames(rt)
		hr <- sprintf("%.3f",rt$"HR")
		hrLow  <- sprintf("%.3f",rt$"HR.95L")
		hrHigh <- sprintf("%.3f",rt$"HR.95H")
		Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
		pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
		pdf(file=forestFile, width = 6.3,height = 4.5)
		n <- nrow(rt)
		nRow <- n+1
		ylim <- c(1,nRow)
		layout(matrix(c(1,2),nc=2),width=c(3,2.5))

		xlim = c(0,3)
		par(mar=c(4,2.5,2,1))
		plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
		text.cex=0.8
		text(0,n:1,gene,adj=0,cex=text.cex)
		text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
		text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
		

		par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
		xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
		plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
		arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
		abline(v=1,col="black",lty=2,lwd=2)
# ⚠ ifelse 的两个分支返回同一个 forestCol, 等价于无条件取该颜色 ——
#   与 44.Prognosis_Related_Gene_Selection 里「HR>1 红 / <1 绿」的做法不同。
#   未自动修改: 若本意就是「每张森林图用单一颜色」, 现状是对的。
		boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
		points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
		axis(1)
		dev.off()
}


bioForest(coxFile="unCox.txt",forestFile="unForest.pdf",forestCol="green")
bioForest(coxFile="muCox.txt",forestFile="muForest.pdf",forestCol="red")
