# 把 12.Expression_and_Immu_cell_activation 的相关性结果画成棒棒糖图
#
# 输入: cor.result.txt(来自 12 章; 需含 Cell / cor / pvalue 三列)
# 输出: corelation.pdf
#
# 本文件不做任何统计, 只负责呈现。

setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/15.Immu_Cell")

# this analysisi is based on the Chapter: 12.Expression_and_Immu_cell_activation

inputFile="cor.result.txt"       #?????ļ?



data = read.table(inputFile, header=T, sep="\t", check.names=F)
data

# 颜色按 pvalue 分档: gold 表示 p>0.8(最不显著), darkgreen 表示 p<=0.2(最显著)。
# 与文件末尾颜色条的刻度(1 → 0)方向一致。
p.col = c('gold','pink','orange','LimeGreen','darkgreen')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                ifelse(x>0.2,p.col[4], p.col[5])
                )))
  return(color)
}

#????????ԲȦ??С?ĺ???
# 圆圈大小按 |cor| 分档: 越大表示相关性越强
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
              ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)
}

#????pvalue????ԲȦ????ɫ
points.color = fcolor(x=data$pvalue,p.col=p.col)
data$points.color = points.color

#????????ϵ??????ԲȦ?Ĵ?С
points.cex = fcex(x=data$cor)
data$points.cex = points.cex
# 按相关系数排序, 使图从下到上单调 —— 纯视觉处理, 不影响统计
data=data[order(data$cor),]

########????ͼ??########
xlim = ceiling(max(abs(data$cor))*10)/10         #x?᷶Χ

pdf(file="corelation.pdf", width=9, height=7)      #????ͼ??
# layout 把画布分成 主图 / 大小图例 / 颜色图例 三块; 数字是写死的版面比例
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(data)+0.5),xlab="Correlation Coefficient",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(data),col="white",lty=1,lwd=2)
#????ͼ?ε??߶?
segments(x0=data$cor,y0=1:nrow(data),x1=0,y1=1:nrow(data),lwd=4)
#????ͼ?ε?ԲȦ
points(x=data$cor,y = 1:nrow(data),col = data$points.color,pch=16,cex=data$points.cex)
#չʾ????ϸ????????
text(par('usr')[1],1:nrow(data),data$Cell,adj=1,xpd=T,cex=1.5)
#չʾpvalue
pvalue.text=ifelse(data$pvalue<0.001,'<0.001',sprintf("%.03f",data$pvalue))
redcutoff_cor=0
redcutoff_pvalue=0.05
# ⚠ redcutoff_cor=0 意味着只要 p<0.05 就标红, 相关系数大小不参与判断
text(par('usr')[2],1:nrow(data),pvalue.text,adj=0,xpd=T,col=ifelse(abs(data$cor)>redcutoff_cor & data$pvalue<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F)

#????ԲȦ??С??ͼ??
par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")

#????ԲȦ??ɫ??ͼ??
par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off()

