library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
library(future.apply)
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
data=data[rowMeans(data)>0.5,]


group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=t(avereps(data))


geneExp=as.data.frame(t(data[gene,,drop=F]))
geneExp$Type=ifelse(geneExp[,gene]>median(geneExp[,gene]), "High", "Low")



for(drug in allDrugs){
  #预测药物敏感性
  possibleError=tryCatch(
    {senstivity=pRRopheticPredict(data, drug, selection=1, dataset = "cgp2016")},
    error=function(e) e)
  if(inherits(possibleError, "error")){next}
  senstivity=senstivity[senstivity!="NaN"]
  senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
  
  #合并敏感性与表达矩阵
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
  outTab=data.frame()
  y=as.numeric(rt$senstivity)
  corT=cor.test(x, y, method = 'pearson')
  cor=corT$estimate
  pvalue=corT$p.value
  outTab=rbind(outTab, cbind(Query=gene, Gene=drug, cor, pvalue))
  #保存满足条件的基因
  #可视化
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
