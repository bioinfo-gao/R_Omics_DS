###Video source: http://study.163.com/u/biowolf
######生信商城：http://www.biowolf.cn/shop/
######速科生物: http://www.biowolf.cn/
######作者QQ：2749657388

#install.packages("pheatmap")

setwd("C:\\Users\\lexb4\\Desktop\\TCGAclinical\\06.pheatmap")      #设置工作目录
rt=read.table("diffmRNAExp.txt",sep="\t",header=T,row.names=1,check.names=F)
rt=log2(rt+1)
rt[rt>15]=15

library(pheatmap)
Type=c(rep("Primary",49),rep("Node",123))    #修改正常和癌症样品数目
names(Type)=colnames(rt)
Type=as.data.frame(Type)

tiff(file="heatmap.tiff",
       width = 45,            #图片的宽度
       height =50,            #图片的高度
       units ="cm",
       compression="lzw",
       bg="white",
       res=300)
pheatmap(rt, annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F,
#         scale="row",
         fontsize_row=5,
         fontsize_col=5)
dev.off()

###Video source: http://study.163.com/u/biowolf
######生信商城：http://www.biowolf.cn/shop/
######速科生物: http://www.biowolf.cn/
######作者QQ：2749657388
