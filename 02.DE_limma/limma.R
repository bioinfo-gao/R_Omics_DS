#biocLite("limma")

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma")
# install.packages(c( "gplots", "beepr", "tidyverse" ))
# install.packages("here")
#library(here) # here 包来构建跨平台的路径，它会自动处理分隔符问题。

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

getwd()
setwd("C:\\Users\\zhen-\\Code\\R_code\\R_For_DS_Omics\\02.差异分析//")
    

rt=read.table("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/01.New_TCGA/combined_RNAseq_FPKM.txt",sep="\t",header=T,check.names=F)

dim(rt)
head(rt)
rt[ 1:3, 1:5 ]

rt=as.matrix(rt)

rownames(rt)=rt[,1]

exp=rt[ ,3:ncol(rt)] # exp=rt[ ,2:ncol(rt)]

exp[1:2, 1:5]

dimnames=list(rownames(exp), colnames(exp))

dimnames # gene names & sample names
?nrow

data=matrix( as.numeric( as.matrix(exp) ),  nrow=nrow(exp),   dimnames = dimnames) # nrow=nrow(exp) <<== is NEEDED !!

data[1:3, 1:3]

# Condense a microarray data object so that values for within-array replicate probes 
# are replaced with their average.
#?avereps()
# Not needed at all here for NGS
# data=avereps(data) 

data=data[rowMeans(data)>0,]

data2=as.data.frame(data)

data2[1:3, 1:3]
colnames(data2)

#??01A??11A???飬??????ǰ?棬?????ź???

exp_data_T = data2 %>% dplyr::select(str_which( colnames(.),  "-01A")) #%>% can NOT be replaced !! 

nT = ncol(exp_data_T) 

nT

exp_data_N = data2 %>% dplyr::select(str_which(colnames(.), "-11A")) # %>% can NOT be replaced 



nN = ncol(exp_data_N) 

nN

rt= cbind(exp_data_N, exp_data_T)


dim(rt)

rt[1:2,1:5]

getwd()
#write.table(rt,file = "groupout.txt",sep="\t",quote=F)
#rt=read.table("groupout.txt",sep="\t",header=T,check.names=F,row.names = 1)

#normalize????Ϊ?????ݵ?У??
pdf(file="rawBox.pdf",  width = 8000) # pixes 


boxplot(rt, col = "blue", xaxt = "n", outline = F) 

dev.off()



rt=normalizeBetweenArrays(rt)

pdf(file="normalBox.pdf")

#??????ͼ
boxplot(rt,col = "red",xaxt = "n",outline = F)
dev.off()

group1=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
group1=sapply(strsplit(group1,""), "[", 1)
group1=gsub("2", "1", group1)

conNum=length(group1[group1==1])       #????????Ʒ??Ŀ
treatNum=length(group1[group1==0])     #????????Ʒ??Ŀ

#differential????????
class <- c(rep("con",conNum),rep("treat",treatNum))  
design <- model.matrix(~factor(class)+0)
colnames(design) <- c("con","treat")

# ==============================================================
df.fit <- lmFit(rt,design)
df.matrix<- makeContrasts(con - treat,levels=design)
fit<- contrasts.fit(df.fit,df.matrix)

#??Ҷ˹????
fit2 <- eBayes(fit)

#????????
allDiff=topTable(fit2,adjust='fdr',n=Inf) 

#д??????
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)


#?页?????两?????希?pvalue小??0.05??写??????
diff
Lab <- allDiff[with(allDiff, ((logFC> 1 | logFC< (-1)) & adj.P.Val < 0.05 )), ]
write.table(diffLab,file="diffExp.xls",sep="\t",quote=F)

#????????????ˮƽ?????ڹ?????
diffExpLevel <- rt[rownames(diffLab),]
write.table(diffExpLevel,file="diffExpLevel.xls",sep="\t",quote=F)

#???ӻ?
gene="THBS2" 
data=t(rt[gene,,drop=F])
Type=c(rep(1,conNum), rep(2,treatNum))
exp=cbind(data, Type)
exp=as.data.frame(exp)
colnames(exp)=c("gene", "Type")
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
exp$gene=log2(exp$gene+1)
group=levels(factor(exp$Type))
exp$Type=factor(exp$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
boxplot=ggboxplot(exp, x="Type", y="gene", color="Type",
                  xlab="",
                  ylab=paste0(gene, " expression"),
                  legend.title="Type",
                  palette = c("blue","red"),
                  add = "jitter")+ 
  stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
#???




?ͼƬ
pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
print(boxplot)
dev.off()

