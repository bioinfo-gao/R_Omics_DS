# 目标基因表达 与 药物敏感性(pRRophetic 预测的 IC50) 的相关性
#
# 输入: combined_RNAseq_counts.txt; 药物响应模型来自 pRRophetic 自带的 CGP2016 数据
# 输出: 仅对达到阈值的药物各出一张 <drug>.pdf 散点图
#       ⚠ 没有输出汇总表, 见第 63 行说明

library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
library(future.apply)
# pRRopheticPredict 内部有随机性, 顶部固定种子是必要的
set.seed(12345)
corFilter=0.3
pFilter=0.05
gene="ARL15"          
expFile="combined_RNAseq_counts.txt"    


data(cgp2016ExprRma)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
allDrugs=unique(drugData2016$Drug.name)


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
# ⚠ 这里仍是【原始 counts】。pRRophetic 的训练集(CGP)是 RMA 芯片数据,
#   预期输入是已归一化并取过对数的表达量。直接喂 counts 会让岭回归外推到
#   完全不同的数值尺度上, 预测出的 IC50 不具可解释性。
#   同目录性质的 52.oncoPredict 脚本就做了 log2(x+1), 两者口径并不一致。
#   未自动修改: 换尺度会改变全部结果, 需你确认。
data=data[rowMeans(data)>0.5,]


group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
# 只保留肿瘤样本, 再把 barcode 截到患者层级
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=t(avereps(data))


geneExp=as.data.frame(t(data[gene,,drop=F]))
geneExp$Type=ifelse(geneExp[,gene]>median(geneExp[,gene]), "High", "Low")



for(drug in allDrugs){
  #预测药物敏感性
# 有些药物在 CGP 里拟合不出模型, 这里用 tryCatch 跳过
  possibleError=tryCatch(
    {senstivity=pRRopheticPredict(data, drug, selection=1, dataset = "cgp2016")},
    error=function(e) e)
  if(inherits(possibleError, "error")){next}
# ⚠ senstivity!="NaN" 是把数值与字符串比较, 靠隐式类型转换才碰巧过滤掉了 NaN。
#   语义明确的写法是 senstivity[!is.nan(senstivity)]。
  senstivity=senstivity[senstivity!="NaN"]
# 把最高 1% 的 IC50 截尾(winsorize), 避免极端值主导相关系数
  senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
  
  #合并敏感性与表达矩阵
# ✅ 这里有按样本名取交集 —— 与 52.oncoPredict 里按位置配对的做法不同, 本文件是对的
  sameSample=intersect(row.names(geneExp), names(senstivity))
  geneExp1=geneExp[sameSample, "Type",drop=F]
  geneExp2=geneExp[sameSample,gene,drop=F]
  senstivity=senstivity[sameSample]
  rt=cbind(geneExp1, senstivity)
  rt$Type=factor(rt$Type, levels=c("Low", "High"))
  type=levels(factor(rt[,"Type"]))
  comp=combn(type, 2)
  
  #提取目标基因表达量
  x=as.numeric(geneExp2[,gene])
# ⚠ outTab 在【循环内部】被重新初始化, 每轮只保留当前药物一行,
#   而且循环结束后从未被 write.table 写出。
#   结果是本脚本只产出图, 不产出任何汇总表。若需要表, 应把这一行移到循环之前。
  outTab=data.frame()
  y=as.numeric(rt$senstivity)
  corT=cor.test(x, y, method = 'pearson')
  cor=corT$estimate
  pvalue=corT$p.value
  outTab=rbind(outTab, cbind(Query=gene, Gene=drug, cor, pvalue))
  #保存满足条件的基因
  #可视化
# ⚠ 逐药物独立检验, 未做多重检验校正; 药物有上百个, 假阳性会很可观
  if((abs(cor)>corFilter) & (pvalue<pFilter)){
    df1=as.data.frame(cbind(x,y))
    p1=ggplot(df1, aes(x, y)) + 
      xlab(paste0(gene, " expression"))+ ylab(paste0(drug, "  drug sensitivity (IC50)"))+
      geom_point()+ geom_smooth(method="lm", formula=y~x) + theme_bw()+
      stat_cor(method = 'pearson', aes(x =x, y =y))
    pdf(file=paste0(drug, ".pdf"), width=5, height=4.6)
    print(p1)
    dev.off()
  }
}
