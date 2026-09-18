###Video source: http://study.163.com/u/biowolf
######生信商城：http://www.biowolf.cn/shop/
######速科生物: http://www.biowolf.cn/
######作者QQ：2749657388

inputFile="clinicalExp.txt"                                     #输入文件
setwd("C:\\Users\\lexb4\\Desktop\\TCGAclinical\\08.ks")         #工作目录
geneName="MPPED1"
clinical="stage"

rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=rt[,c("id",clinical,geneName)]
colnames(rt)=c("id","clinical","gene")

xlabel=vector()
tab1=table(rt[,"clinical"])
labelNum=length(tab1)
dotCol=c(2,3)
if(labelNum==3){
	dotCol=c(2,3,4)
}
if(labelNum==4){
	dotCol=c(2,3,4,5)
}
if(labelNum>4){
	dotCol=rainbow(labelNum)
}
for(i in 1:labelNum){
  xlabel=c(xlabel,names(tab1[i]))
}

ksTest<-kruskal.test(gene ~ clinical, data = rt)
wilcoxP=ksTest$p.value
pvalue=signif(wilcoxP,4)
pval=round(pvalue,3)

b = boxplot(gene ~ clinical, data = rt,outline = FALSE, plot=F) 
yMin=min(b$stats)
yMax = max(b$stats/5+b$stats)
ySeg = max(b$stats/10+b$stats)
ySeg2 = max(b$stats/12+b$stats)
n = ncol(b$stats)

tiffFile=paste(geneName,".tiff",sep="")
tiff(file=tiffFile,width = 26,height = 16,
     units ="cm",compression="lzw",bg="white",res=300)
par(mar = c(4,7,3,3))
boxplot(gene ~ clinical, data = rt,names=xlabel,
     ylab = paste(geneName," expression",sep=""),col=dotCol,
     cex.main=1.6, cex.lab=1.4, cex.axis=1.3,ylim=c(yMin,yMax),outline = FALSE)
segments(1,ySeg, n,ySeg);
segments(1,ySeg, 1,ySeg2)
segments(n,ySeg, n,ySeg2)
text((1+n)/2,ySeg,labels=paste("p=",pval,sep=""),cex=1.5,pos=3)
dev.off()

###Video source: http://study.163.com/u/biowolf
######生信商城：http://www.biowolf.cn/shop/
######速科生物: http://www.biowolf.cn/
######作者QQ：2749657388
