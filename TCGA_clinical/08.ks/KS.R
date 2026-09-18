# 单基因表达在不同临床分期(或其它分类变量)间的差异 —— Kruskal-Wallis 检验 + 箱线图
#
# 输入: clinicalExp.txt (需含 id 列、指定的临床列、以及目标基因列)
# 输出: <geneName>.tiff
#
# 注: 本文件来自第三方教学材料(biowolf), 非本项目自写。

###Video source: http://study.163.com/u/biowolf
######生信商城：http://www.biowolf.cn/shop/
######速科生物: http://www.biowolf.cn/
######作者QQ：2749657388

inputFile="clinicalExp.txt"                                     #输入文件
# ⚠ setwd 指向教程作者的机器路径, 运行前需改
setwd("C:\\Users\\lexb4\\Desktop\\TCGAclinical\\08.ks")         #工作目录
geneName="MPPED1"
clinical="stage"

rt=read.table(inputFile,sep="\t",header=T,check.names=F)
# ⚠ 按列名取三列; clinicalExp.txt 里必须确实存在 id / <clinical> / <geneName> 这三列,
#   否则这一行直接报 undefined columns。
rt=rt[,c("id",clinical,geneName)]
colnames(rt)=c("id","clinical","gene")

# 下面这段只是根据分组数挑配色与刻度标签, 与统计无关
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

# ⚠ 两点:
#   1) 这里用的是 kruskal.test(多组), 但下面把结果存进了名为 wilcoxP 的变量 —— 命名会误导;
#   2) 未剔除缺失值或 "unknow" 这类占位分组。本仓库的 06_clinical / 70.Clinical_* 都做了
#      data[data$clinical!="unknow",] 的过滤, 这里没有, 占位值会被当成一个真实的组参与检验。
ksTest<-kruskal.test(gene ~ clinical, data = rt)
wilcoxP=ksTest$p.value
pvalue=signif(wilcoxP,4)
pval=round(pvalue,3)

# 先跑一次不绘图的 boxplot 取分位数, 用来给显著性横线留出纵向空间
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
# ⚠ outline=FALSE 只是不【画】离群点, 它们仍参与上面的检验 —— 图与统计的样本并不一致, 读图时留意
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
