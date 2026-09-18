# 目标基因表达 与 各免疫细胞浸润比例 的相关性 (Spearman)
#
# 输入: geneExp.txt(需含 Type 列, 第 1 列为目标基因), CIBERSORT-Results.txt
# 输出: cor.result.txt(供 15.Immu_Cell_Activation_Exp_Correlation 画图),
#       以及每个显著细胞一张 cor.<cell>.pdf
#
# 注: 本章的后续分析在 15.Immu_Cell_Activation_Exp_Correlation_Ch12。

setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/12.Expression_and_Immu")
# The chapter is followed by Chapter 15 <<==============
library(limma)
library(reshape2)
library(ggpubr)
library(vioplot)
library(ggExtra)

expFile="geneExp.txt"              #?????????ļ?
immFile="CIBERSORT-Results.txt"    #????ϸ???????Ľ????ļ?
pFilter=0.05            #????ϸ???????????Ĺ???????


#??ȡ?????????ļ?
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
# ⚠ 硬编码取第 1 列作为目标基因
gene=colnames(rt)[1]

#ɾ????????Ʒ
tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

#????目??????????量????品???蟹???
data=as.data.frame(data)
data$gene=ifelse(data[,gene]>median(data[,gene]), "High", "Low")

#??ȡ????ϸ???????ļ??????????ݽ???????
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
# CIBERSORT 的 P-value 衡量的是该样本去卷积整体是否可信, 不是单个细胞类型的显著性
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

#ɾ????????Ʒ
group=sapply(strsplit(row.names(immune),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
immune=immune[group==0,]
row.names(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(immune))
immune=avereps(immune)

#???ݺϲ?
sameSample=intersect(row.names(immune), row.names(data))
rt=cbind(immune[sameSample,,drop=F], data[sameSample,,drop=F])


# ⚠ 下面这三行(melt 出来的 data)在本文件中【从未被使用】—— 是从 10.Immu_cell 抄来的残留;
#   真正参与循环的是上面的 rt。
#??????ת????ggplot2?????ļ?
data=rt[,-(ncol(rt)-1)]
data=melt(data,id.vars=c("gene"))
colnames(data)=c("gene", "Immune", "Expression")
##########??????????ɢ??ͼ##########
outTab=data.frame()
# 排除最后两列(目标基因表达值与 High/Low 分组), 其余都是免疫细胞比例列
for(i in colnames(rt)[1:(ncol(rt)-2)]){
  x=as.numeric(rt[,gene])
  y=as.numeric(rt[,i])
# ⚠ 这是为了绕开「方差为 0 时 cor.test 报错」而人为改动一个观测值。
#   它会让该细胞类型的相关系数变成一个无意义的数, 而不是如实报告「无法计算」。
  if(sd(y)==0){y[1]=0.00001}
# ⚠ Spearman 在有大量并列值(比例常含很多 0)时会告警且 p 值用近似算法;
#   免疫浸润数据尤其容易触发。
  cor=cor.test(x, y, method="spearman")
  outVector=cbind(Cell=i, cor=cor$estimate, pvalue=cor$p.value)
  outTab=rbind(outTab,outVector)
# ⚠ 未做多重检验校正, 且只有 p<0.05 的才出图
  if(cor$p.value<0.05){
    outFile=paste0("cor.", i, ".pdf")
    df1=as.data.frame(cbind(x,y))
    p1=ggplot(df1, aes(x, y)) + 
      xlab(paste0(gene, " expression")) + ylab(i)+
      geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
      stat_cor(method = 'spearman', aes(x =x, y =y))
    p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
    #??????ͼ??
    pdf(file=outFile, width=5.2, height=5)
    print(p2)
    dev.off()
  }
}
#?????????ԵĽ????ļ?
write.table(outTab,file="cor.result.txt",sep="\t",row.names=F,quote=F)

