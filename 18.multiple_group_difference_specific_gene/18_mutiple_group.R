# 多分组(A/B/C/D 四个文件)合并后, 比较单个基因在各组间的表达差异
#
# 输入: A.txt / B.txt / C.txt / D.txt, 均为 行=基因 列=样本, 第 1 列是基因名
# 输出: all.txt(合并矩阵), groups.txt(样本-分组对照), TypesCor_group.pdf(箱线图)
#
# ⚠⚠ 本文件此前被一次有损编码转换破坏: 共 6 行代码被并进注释, 导致
#   groups / samSample / gene / group / boxplot 全部未定义, 且图不落盘。
#   六处均已还原, 每处在原地注明判据。

setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/18.multiple_group_difference")
library(limma)
library(ggpubr)

rt1=read.table("A.txt",sep="\t",header=T,check.names=F)
rt2=read.table("B.txt",sep="\t",header=T,check.names=F)
rt3=read.table("C.txt",sep="\t",header=T,check.names=F)
rt4=read.table("D.txt",sep="\t",header=T,check.names=F)
# rt1=as.matrix(rt1)
# rt2=as.matrix(rt2)
# rt3=as.matrix(rt3)
# rt4=as.matrix(rt4)
rownames(rt1)=rt1[,1]
rownames(rt2)=rt2[,1]
rownames(rt3)=rt3[,1]
rownames(rt4)=rt4[,1]
rt1=rt1[,2:ncol(rt1)]
rt2=rt2[,2:ncol(rt2)]
rt3=rt3[,2:ncol(rt3)]
rt4=rt4[,2:ncol(rt4)]

# 四个文件按列拼接; ⚠ 依赖四者【行名(基因)完全一致且顺序相同】, cbind 不做任何对齐检查
bindrt=cbind(rt1,rt2,rt3,rt4)
dimnames=list(rownames(bindrt),colnames(bindrt))
data=matrix(as.numeric(as.matrix(bindrt)),nrow=nrow(bindrt),dimnames=dimnames)
write.table(data, file="all.txt", sep="\t", quote=F, row.names=T)

# 分组标签取自列名的第 1 个字符 —— 即样本名首字母就是组别(A/B/C/D)
# 已还原: 原为 groups=sapply(...), 变量名被编码事故截成了 ps; 判据是下一行使用 groups。
groups=sapply(strsplit(colnames(data),""), "[", 1)
Type=as.matrix(groups)
colnames(Type)="group"
ids=as.matrix(colnames(as.matrix(data)))
colnames(ids)="id"
groupbind=cbind(ids, groups)
groupbind=as.matrix(groupbind)
rownames(groupbind)=groupbind[,1]
groupbind=as.data.frame(groupbind[,2:ncol(groupbind)])
colnames(groupbind)="group"
write.table(groupbind, file="groups.txt", sep="\t", quote=F, row.names=T)

# 注: data[drop=F,,samSample] 是把 drop 作为命名参数先匹配, 余下两个位置参数
#     分别落到 i(空) 和 j, 等价于 data[, samSample, drop=F]。写法古怪但有效。
# 已还原: 原为 samSample=intersect(...), 整行被并进了注释; 判据是下一行使用 samSample。
samSample=intersect(colnames(data), row.names(groupbind))
data=t(data[drop=F,,samSample])
cli=groupbind[samSample,,drop=F]
rt=cbind(data, cli)

# 已还原: 原为 gene="TSPAN6"(目标基因), 整行被并进注释; 判据是下一行使用 gene。
gene="TSPAN6"
# 剔除 "unknow"(原拼写)分组
Types="group"
data=rt[c(gene, Types)]
colnames(data)=c(gene, "Types")
data=data[(data[,"Types"]!="unknow"),]
# 已还原: 原为 group=levels(factor(data$Types)); 判据是下面两行使用 group。
group=levels(factor(data$Types))
data$Types=factor(data$Types, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# ⚠ stat_compare_means 默认做的是两两比较, 未对多组比较做校正;
#   组数一多, 星号里会出现假阳性。
# 已还原: 原为 boxplot=ggboxplot(...); 判据是文件末尾 print(boxplot)。
boxplot=ggboxplot(data, x="Types", y=gene, fill="Types",
                  xlab=Types,
                  ylab=paste(gene, " expression"),
                  legend.title=Types)+ 
  stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
# 已还原: 原为 pdf(file=...); 被吞进注释后图不落盘, 且下面的 dev.off() 关的是默认设备。
pdf(file=paste0("TypesCor_", Types, ".pdf"), width=5.5, height=5)
print(boxplot)
dev.off()
