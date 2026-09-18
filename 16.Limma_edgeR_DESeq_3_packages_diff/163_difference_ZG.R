# 163_difference.R 的修整版: 同一套数据用 limma / DESeq2 / edgeR 三法比较
#
# 输入/输出与 163_difference.R 相同
#
# 本文件是两份副本中较完整的一份: install.packages 已注释、韦恩图标题是英文,
# 且 data=colnames(countData) 与 edgeR 段的 read.table 都完好 ——
# 这三处正是 163_difference.R 中被编码事故破坏、后来据本文件还原的地方。
#
# ⚠ 与另一份相同的问题: edgeR 段阈值原本与前两段不一致, 现已统一。

setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/16.Limma_edgeR_DESeq_3_packages_diff")
#install.packages("VennDiagram")
library(edgeR)
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
library("DESeq2")
library(VennDiagram)

padj = 0.05
foldChange= 2


##### (1) limma DE #####  <<<===============================================

rt1=read.table("../01.New_TCGA//combined_RNAseq_FPKM.txt",sep="\t",header=T,check.names=F)
dim(rt1)
rt1[1:2, 1:5]

rt1=as.matrix(rt1)
class(rt1)

# SEE THE 41 Utilities ==> gene name

# # 设置行名为唯一基因名，假设基因名在第一列。
# # 注意：如果文件第一列是基因名，直接用 row.names=1 更好，但为了保持原代码逻辑，我们使用 make.unique。
# # 此处假定 rt[,1] 是基因名
# rownames(rt) = make.unique(rt[,2]) 
# 
# # 提取表达数据。原代码注释提到 "not 2"，通常基因名是第1列。
# # 故从第2列（或第3列，取决于文件结构）开始取表达值。这里沿用原代码的 rt[,3:ncol(rt)]。
# # 如果第2列是基因 ID，则从第3列开始；如果第2列就是第一个样本，则应该从 rt[,2:ncol(rt)]。

rownames(rt1)=make.unique(rt1[,1])

exp1=rt1[,3:ncol(rt1)] # exp1=rt1[,2:ncol(rt1)]
exp1[1:2, 1:5]

dimnames=list(rownames(exp1),colnames(exp1))
data1=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames)
data1=avereps(data1)
data1=data1[rowMeans(data1)>0,]
data1=as.data.frame(data1)

dim(data1)
data1[1:2, 1:5]


# ===========================<<<<<<<<<<<<<<<

exp1_data_T = data1%>% dplyr::select(str_which(colnames(.), "-01A")) # 
nT = ncol(exp1_data_T) 
exp1_data_N = data1%>% dplyr::select(str_which(colnames(.), "-11A"))
nN = ncol(exp1_data_N) 
rt1= cbind(exp1_data_N, exp1_data_T)

rt1=normalizeBetweenArrays(rt1)
group1=sapply(strsplit(colnames(rt1),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1) # No seperator means split string to character
group1=gsub("2", "1", group1)

conNum1=length(group1[group1==1])       #????????Ʒ??Ŀ
treatNum1=length(group1[group1==0])     #????????Ʒ??Ŀ

#differential Analysis
class <- c(rep("con",conNum1),rep("treat",treatNum1))  
design <- model.matrix(~factor(class)+0)
colnames(design) <- c("con","treat")

df.fit <- lmFit(rt1,design)
df.matrix<- makeContrasts(con - treat,levels=design)
fit<- contrasts.fit(df.fit,df.matrix)
#??Ҷ˹????
fit2 <- eBayes(fit)
#????????
allDEG1 = topTable(fit2,coef=1,n=Inf,adjust="BH") 
allDEG1 = na.omit(allDEG1)



padj = 0.05
foldChange= 2
diff_signif1 = allDEG1[(allDEG1$adj.P.Val < padj & 
                          (allDEG1$logFC>foldChange | allDEG1$logFC<(-foldChange))),]
diff_signif1 = diff_signif1[order(diff_signif1$logFC),]
save(diff_signif1, file = 'limma_diff1.Rdata')



##### (2) DESeq DE ##### <<<===============================================
rt <- read.table( "../01.New_TCGA/combined_RNAseq_counts.txt",header=T,sep="\t",comment.char="",check.names=F)
rt=as.matrix(rt)

rt[1:2, 1:5]

rownames(rt)= make.unique(rt[,1]) 

exp=rt[,3:ncol(rt)] # 3

dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]
data2=as.data.frame(data)
data2[1:2, 1:5]

exp_data_T = data2 %>% dplyr::select(str_which(colnames(.), "-01A$")) # ƥ????????????ʾд??
nT = ncol(exp_data_T) 
exp_data_N = data2 %>% dplyr::select(ends_with("-11A"))
nN = ncol(exp_data_N) 
data= cbind(exp_data_N, exp_data_T)


group1=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)
conNum=length(group1[group1==1])       #????????Ʒ??
treatNum=length(group1[group1==0])     #????????Ʒ??

count <- floor(data) # EDSeg data is continuous 

countData <- count[apply(count, 1, sum) > 0 , ]

data=colnames(countData)
Type=c(rep(1,conNum), rep(2,treatNum))
exp=cbind(data, Type)
exp=as.data.frame(exp)
colnames(exp)=c("id", "Type")
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
exp=as.matrix(exp)
rownames(exp)=exp[,1]
exp=exp[,2:ncol(exp)]
exp=as.data.frame(exp)
colnames(exp)=c("condition")
write.table(exp, file="group.txt",sep="\t",quote=F)
colData <- read.table( "group.txt",header=T,sep="\t",row.names=1,comment.char="",check.names=F)


dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "Normal")
dds$condition
summary(dds)
dds

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)  # 15 min 
# gene-wise dispersion estimates
# gene-wise dispersion estimates 
# mean-dispersion relationship
# final dispersion estimates

dds <- nbinomWaldTest(dds)
dds <- DESeq(dds)              # another 18 min
allDEG2 <- as.data.frame(results(dds))

padj = 0.05
foldChange= 2
diff_signif2 = allDEG2[(allDEG2$padj < padj & (allDEG2$log2FoldChange > foldChange | allDEG2$log2FoldChange < (-foldChange))),]
diff_signif2 = diff_signif2[order(diff_signif2$log2FoldChange),]
save(diff_signif2, file = 'DESeq2_diff2.Rdata')


##### (3) degeR DE #####  <<<=============================================== HERE 
rt=read.table("../01.New_TCGA/combined_RNAseq_counts.txt",sep="\t",header=T,check.names=F) 
#rt=as.matrix(rt)  
rt[1:2, 1:5]

rownames(rt)=make.unique(rt[,1])


exp=rt[,3:ncol(rt)] # changed from 2 
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]
data2=as.data.frame(data)
exp_data_T = data2%>% dplyr::select(str_which(colnames(.), "-01A$")) 
nT = ncol(exp_data_T) 
exp_data_N = data2%>% dplyr::select(ends_with("-11A"))
nN = ncol(exp_data_N) 
data= cbind(exp_data_N, exp_data_T)
group1=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)
conNum=length(group1[group1==1])       
treatNum=length(group1[group1==0])     

group=c(rep("normal",conNum),rep("tumor",treatNum)) 

design <- model.matrix(~group)

y <- DGEList(counts=data,group=group)#?????б?
y <- calcNormFactors(y)#?????????ڱ?׼??????
y <- estimateCommonDisp(y)#??????ͨ????ɢ
y <- estimateTagwiseDisp(y)#??????????miRNA??Χ?ڵ???ɢ
et <- exactTest(y,pair = c("normal","tumor"))#???о?ȷ??

topTags(et) 

ordered_tags <- topTags(et, n=100000) 

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]

# ⚠ 已改: 原为 foldChange=1, 而上面 limma 与 DESeq2 两段用的是 2。
#   三种方法必须用同一条尺子, 否则末尾的韦恩图比较的是不同阈值下的结果。
padj = 0.05
foldChange= 2
diff_signif = allDiff[(allDiff$FDR < padj & (allDiff$logFC>foldChange | allDiff$logFC<(-foldChange))),]
diff_signif = diff_signif[order(diff_signif$logFC),]
save(diff_signif, file = 'edger_diff.Rdata')


#####  ==== Difference among the 3 packages ====  #####
edgeR = rownames(diff_signif)
dim(diff_signif)
limma = rownames(diff_signif1)
dim(diff_signif1)
DESeq2 = rownames(diff_signif2)
dim(diff_signif2)

venn.diagram(
  x = list(
    'edgeR' = edgeR,
    'limma' = limma,
    'DESeq2' = DESeq2
  ),
  filename = 'VN.png',
  col = "black",
  fill = c("dodgerblue", "goldenrod1", "darkorange1"),
  alpha = 0.5,
  cex = 0.8,
  cat.col = 'black',
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.05,
  main = "The DE gene detected by Limma edgeR and DESeq2",
  main.cex = 1.2
)

#####????????#####
edgeR=as.data.frame(edgeR)
limma=as.data.frame(limma)
DESeq2=as.data.frame(DESeq2)
sameSample=intersect(edgeR$edgeR, limma$limma)
sameSample=as.data.frame(sameSample)
sameSample=intersect(DESeq2$DESeq2,sameSample$sameSample)

sameSample1=as.data.frame(sameSample)
write.table(sameSample1, file="merge_genes.xls",sep="\t",quote=F)

data1=exp[sameSample,,drop=F]
write.table(data1, file="merge_counts.xls",sep="\t",quote=F)

data2=exp1[sameSample,,drop=F]
write.table(data2, file="merge_FPKM.xls",sep="\t",quote=F)



