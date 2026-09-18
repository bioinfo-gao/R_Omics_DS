# TCGA 与 GEO 两套数据各自做差异分析, 取交集基因并画韦恩图
#
# 输入: combined_RNAseq_FPKM.txt (来自 01.New_TCGA)
#       geo_exp.csv   行=基因, 列=样本
#       group.csv     两列: 样本名, 分组标签(N / T)
# 输出: TCGA_limmaTab.xls / TCGA_difflab.xls / geo_deg_all.txt / geo_upanddown.csv
#       TCGAandGEO_venn.pdf / TCGAandGEO_Genes.txt / TCGA_deg_exp.csv / GEO_deg_exp.csv
#
# ⚠ 两套数据的 logFC 正负号定义相反, 详见下方 makeContrasts 与 coef=2 两处注释

library(venn) 
library(VennDiagram)
library(data.table)
library(tidyverse)
library(ggsignif) 
library(RColorBrewer)
library(limma)
library(ggplot2)
library(ggpubr)
library(beepr)
library(gplots)
library(pheatmap)

#####处理TCGA数据#####
rt=read.table("combined_RNAseq_FPKM.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data2=as.data.frame(data)
# TCGA barcode 第 4 段: 01A=原发肿瘤, 11A=癌旁正常
exp_data_T = data2%>% dplyr::select(str_which(colnames(.), "-01A"))
nT = ncol(exp_data_T) 
exp_data_N = data2%>% dplyr::select(str_which(colnames(.), "-11A"))
nN = ncol(exp_data_N) 
# ⚠ 列顺序在此被固定为「先全部正常, 后全部肿瘤」。
#   下面的 class 向量用 rep 按位置连续填充, 不是按列名逐列匹配。
#   这一行的顺序一旦改动, 分组标签会静默反转且不报任何错。
rt= cbind(exp_data_N, exp_data_T)
# 对 FPKM 再做分位数标准化: FPKM 已校正长度与深度, 这一步是为拉齐样本间分布
rt=normalizeBetweenArrays(rt)
# 从 barcode 取样本类型: 第 4 段首字符, 0=肿瘤 1=正常; gsub 把 2(复发)并入 1
# ⚠ group1 此后并未参与分组, 只用于计数 —— 真正决定标签的是上面的列顺序
group1=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)
conNum=length(group1[group1==1])     
treatNum=length(group1[group1==0])     
class <- c(rep("con",conNum),rep("treat",treatNum))  
design <- model.matrix(~factor(class)+0)
colnames(design) <- c("con","treat")
df.fit <- lmFit(rt,design)
# ⚠ 方向: con - treat = 正常 - 肿瘤, 故 logFC > 0 表示在【正常】中更高。
#   这与下面 GEO 部分的方向相反, 见 coef=2 处。
df.matrix<- makeContrasts(con - treat,levels=design)
fit<- contrasts.fit(df.fit,df.matrix)
fit2 <- eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',n=Inf) 
write.table(allDiff,file="TCGA_limmaTab.xls",sep="\t",quote=F)
# 阈值 |logFC|>1 且 adj.P<0.05; 取绝对值, 故与上下调方向无关
TCGA_diffLab <- allDiff[with(allDiff, ((logFC> 1 | logFC< (-1)) & adj.P.Val < 0.05 )), ]
TCGA_genes<-as.data.frame(row.names(TCGA_diffLab))
colnames(TCGA_genes)="genes"
rownames(TCGA_genes)=TCGA_genes$genes
write.table(TCGA_diffLab,file="TCGA_difflab.xls",sep="\t",quote=F)

#####处理GEO数据#####
geo=read.csv("geo_exp.csv",header=T,check.names=F)
geo=as.matrix(geo)
rownames(geo)=geo[,1]
exp_geo=geo[,2:ncol(geo)]
dimnames_geo=list(rownames(exp_geo),colnames(exp_geo))
data_geo=matrix(as.numeric(as.matrix(exp_geo)),nrow=nrow(exp_geo),dimnames=dimnames_geo)
data_geo=avereps(data_geo)
data_geo=data_geo[rowMeans(data_geo)>0,]
data2_geo=as.data.frame(data_geo)
data2_geo=normalizeBetweenArrays(data2_geo)
group_geo=read.csv("group.csv")
group1_geo=as.data.frame(group_geo[,2:ncol(group_geo)])
rownames(group1_geo)=group_geo[,1]
colnames(group1_geo)="groups"
group1_geo_list=factor(group1_geo$groups,levels = c("N","T"))
# ⚠ levels=c("N","T") 令 N 为参照, coef=2 取的是 T vs N,
#   故 GEO 的 logFC > 0 表示在【肿瘤】中更高 —— 与上面 TCGA 正好相反。
#   下游取交集只用 |logFC| 阈值, 方向被忽略, 所以韦恩图是「方向无关」的交集。
#   若要方向一致的共同 DEG, 必须先把其中一套的 logFC 取反再比较。
design=model.matrix(~group1_geo_list)
fit=lmFit(data2_geo,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
write.table(deg, file = "geo_deg_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
k1 = (deg$adj.P.Val < 0.05)&(deg$logFC < -1)
k2 = (deg$adj.P.Val < 0.05)&(deg$logFC > 1)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
write.csv(deg,file="geo_upanddown.csv")
geo_diffLab <- deg[with(deg, ((logFC> 1 | logFC< (-1)) & adj.P.Val < 0.05 )), ]
geo_genes<-as.data.frame(row.names(geo_diffLab))
colnames(geo_genes)="genes"
rownames(geo_genes)=geo_genes$genes

#####提取共有基因绘制韦恩图#####
# 交集只基于基因名, 不区分上调/下调方向(原因见上)
genesList=list()
geneNames_TCGA=as.vector(TCGA_genes[,1])       
uniq_TCGAGene=unique(geneNames_TCGA)       
genesList[["TCGA"]]=uniq_TCGAGene
geneNames_GEO=as.vector(geo_genes[,1])       
uniq_GEOGene=unique(geneNames_GEO)       
genesList[["GEO"]]=uniq_GEOGene
mycol=c("blue","red","yellow","green","orange")
venn.plot <-venn.diagram(genesList,
                         fill=mycol[1:length(genesList)], 
                         filename=NULL, 
                         cat.pos = c(360, 360),
                         scaled =FALSE,#圆圈大小是否按数字大小改变
                         cex = 2,        #1 2 区域内部数字的字体大小 
                         cat.cex = 2,    # 分类名称的字体大小 
                         cat.dist = 0.015,   #分类名称距离边的距离 实际调整 
                         )
pdf("TCGAandGEO_venn.pdf");
grid.draw(venn.plot);
dev.off()
intersectGenes=Reduce(intersect,genesList)
write.table(file="TCGAandGEO_Genes.txt", intersectGenes, sep="\t", quote=F, col.names=F, row.names=F)


#####分别提取共同差异基因的表达矩阵#####
# 此处 rt 已是 normalizeBetweenArrays 的输出(正常列在前, 肿瘤列在后), 非最初读入的 rt
TCGA_deg_exp=rt[intersectGenes,,drop=F]
GEO_deg_exp=data2_geo[intersectGenes,,drop=F]
write.csv(file="TCGA_deg_exp.csv", TCGA_deg_exp)
write.csv(file="GEO_deg_exp.csv", GEO_deg_exp)

