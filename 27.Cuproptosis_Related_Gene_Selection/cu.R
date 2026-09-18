# 筛选与铜死亡(cuproptosis)相关基因: 先做肿瘤/正常差异, 再与铜死亡基因求相关
#
# 输入: combined_RNAseq_FPKM.txt, gene.txt(铜死亡基因列表, 单列无表头)
# 输出: cuproptosis_gene_exp.txt(铜死亡基因表达), Result.txt(显著相关的基因对),
#       cuproptosisExp.txt(相关 RNA 的表达矩阵)
#
# ⚠⚠ 第 33 行的分组标签与列顺序不匹配, 会让第 57 行的差异检验比错组, 见该处说明。


library(limma)            #???ð?
expFile="combined_RNAseq_FPKM.txt"      #?????????ļ?
geneFile="gene.txt"       #?????б??ļ?


#??ȡ?????ļ??????????ݽ??д???
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#??取?????斜??募?,??取铜???????鼗????谋???量
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
# 只取表达矩阵中确实存在的铜死亡基因
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

#????????
outTab=rbind(ID=colnames(geneExp),geneExp)
# 注: 原文件名含不可还原的乱码字符, 已改为 ASCII 名; 第 36 行的读取同步改过
write.table(outTab, file="cuproptosis_gene_exp.txt", sep="\t", quote=F, col.names=F)

#####??????????
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
# RNA 只保留肿瘤样本
RNA=data[,group==0]
conNum=length(group[group==1])       #????????Ʒ??Ŀ
treatNum=length(group[group==0])     #????????Ʒ??Ŀ
# ⚠⚠ sampleType 按「先 conNum 个 1, 再 treatNum 个 2」顺序构造,
#   这隐含假设 data 的列是【正常在前、肿瘤在后】排好的。
#   但本脚本从第 8 行读入后【从未对列做过重排】—— 对比 13.edgeR / 14.DESeq2 /
#   59.GEO_TCGA_Common_DEGs, 那几个脚本都有一步 cbind(exp_data_N, exp_data_T)。
#   因此这里的标签是按文件原始列序贴上去的, 与真实样本类型很可能对不上,
#   下面第 57 行的 Wilcoxon 检验会比错组, 全程不报错。
#   正确做法: 直接用上面已经算好的 group 向量(0=肿瘤 1=正常)构造标签,
#   例如 sampleType = ifelse(group=="1", 1, 2)。
#   未自动修改: 需先确认你的输入文件列序, 以免把本来正确的结果改坏。
sampleType=c(rep(1,conNum), rep(2,treatNum))


rt1=read.table("cuproptosis_gene_exp.txt", header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
cuproptosis=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
cuproptosis=avereps(cuproptosis)
cuproptosis=cuproptosis[rowMeans(cuproptosis)>0,]

#ɾ????????Ʒ
group=sapply(strsplit(colnames(cuproptosis),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
cuproptosis=cuproptosis[,group==0]

#?????Լ???
corFilter=0.3           
pvalueFilter=0.05  
outTab=data.frame()
# ⚠ 双重循环: 外层遍历【全部基因】, 内层遍历铜死亡基因。
#   两万基因 × 十几个铜死亡基因的 cor.test 会很慢, 且未做多重检验校正。
for(i in row.names(RNA)){
# 只检验表达有波动的基因; 注意这里用 RNA(仅肿瘤)判断方差, 下一行却用 data(全部样本)做检验
  if(sd(RNA[i,])>0.1){
# ⚠ 承上: 方差筛选基于 RNA(肿瘤), 差异检验却基于 data(含正常), 两者样本集不一致
    test=wilcox.test(data[i,] ~ sampleType)
    if(test$p.value<0.05){
      for(j in row.names(cuproptosis)){
        x=as.numeric(RNA[i,])
        y=as.numeric(cuproptosis[j,])
        corT=cor.test(x,y)
        cor=corT$estimate
        pvalue=corT$p.value
# 正相关与负相关分开记录, 便于下游按方向解读
        if((cor>corFilter) & (pvalue<pvalueFilter)){
          outTab=rbind(outTab,cbind(Cuproptosis=j,RNA=i,cor,pvalue,Regulation="postive"))
        }
        if((cor< -corFilter) & (pvalue<pvalueFilter)){
          outTab=rbind(outTab,cbind(Cuproptosis=j,RNA=i,cor,pvalue,Regulation="negative"))
        }
      }
    }
  }
}

#?????????ԵĽ???
write.table(file="Result.txt",outTab,sep="\t",quote=F,row.names=F)

#??取铜????????RNA?谋???量
cuproptosisRNA=unique(as.vector(outTab[,"RNA"]))
# ⚠ 这里取的是 data(含正常样本)的表达, 而上面的相关性是在 RNA(仅肿瘤)上算的,
#   导出的表达矩阵与产生它的分析所用样本集并不一致。
cuproptosisRNAexp=data[cuproptosisRNA,]
cuproptosisRNAexp=rbind(ID=colnames(cuproptosisRNAexp), cuproptosisRNAexp)
write.table(cuproptosisRNAexp,file="cuproptosisExp.txt",sep="\t",quote=F,col.names=F)
