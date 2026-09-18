# 单因素 Cox 预后筛选 + 森林图
#
# 输入: exptime.txt —— 列 1=样本id(→rownames) 2=futime 3=fustat 4+=基因表达
# 输出: UniCox.txt(HR/95%CI/p), UniSigExp.txt(通过筛选的基因表达, 供 46/71 建模),
#       unforest.pdf(森林图)


library(survival)
pFilter=0.05                                                       
rt=read.table("exptime.txt",header=T,sep="\t",check.names=F,row.names=1)    

outTab=data.frame()
sigGenes=c("futime","fustat")
# 3:ncol 按位置取基因列, 依赖 futime/fustat 恰在前两列
for(i in colnames(rt[,3:ncol(rt)])){
# 注: 这里把 rt[,i] 直接写进公式, 因此模型里的系数名会是 "rt[, i]" 而不是基因名;
#     下面是用循环变量 i 另行记录基因名的, 所以结果表不受影响。
 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
 coxSummary = summary(cox)
 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
# ⚠ 用的是未校正的原始 p。基因数上千时, p<0.05 会保留大量假阳性,
#   而这一步的输出正是后面 LASSO 建模的候选集 —— 噪声会一路带下去。
#   若要收紧, 可改用 p.adjust(..., method="BH") 后再筛。
 if(coxP<pFilter){
     sigGenes=c(sigGenes,i)
		 outTab=rbind(outTab,
		              cbind(id=i,
		              HR=coxSummary$conf.int[,"exp(coef)"],
		              HR.95L=coxSummary$conf.int[,"lower .95"],
		              HR.95H=coxSummary$conf.int[,"upper .95"],
		              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
		              )
  }
}
write.table(outTab,file="UniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="UniSigExp.txt",sep="\t",row.names=F,quote=F)


######绘制森林图######
#读取输入文件
rt <- read.table("UniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#输出图形
# 以下为手工绘制的森林图: 左半幅放基因名/p/HR 文本, 右半幅放置信区间线段
pdf(file="unforest.pdf", width = 6,height = 4.5)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

#绘制森林图左边的基因信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#绘制森林图
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
# ⚠ x 轴上限只取自 HR 的上下限, 未考虑 HR 本身可能超出该范围(通常不会)
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
# HR>1 标红(危险因素), HR<1 标绿(保护因素)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

