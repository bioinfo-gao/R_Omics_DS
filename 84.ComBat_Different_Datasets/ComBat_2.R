

#引用包
library(limma)
library(sva)
library(FactoMineR)
library(factoextra)
#读取基因表达文件,并对数据进行处理
rt=read.table("combined_RNAseq_FPKM.txt", header=T, sep="\t", check.names=F)
data=as.data.frame(rt)
data0=avereps(data)
#读取GEO数据
GEO=read.csv("geo_exp.csv", header=T, sep=",", check.names=F)
rt=as.matrix(GEO)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data2=as.data.frame(data)
#合并TCGA和GEO数据
sameGene=intersect(rownames(data0),rownames(data2))
TCGA=data0[sameGene,,drop=F]
GEO=data2[sameGene,,drop=F]
##TCGA和GEO组间矫正
TCGA=normalizeBetweenArrays(TCGA)
GEO=normalizeBetweenArrays(GEO)
###组间校正后数据合并
data0=cbind(TCGA,GEO)

#构建批次信息
TCGA_smaple=as.data.frame(colnames(TCGA))
colnames(TCGA_smaple)="group"
TCGA_smaple$batch="TCGA"
GEO_sample=as.data.frame(colnames(GEO))
colnames(GEO_sample)="group"
GEO_sample$batch="GEO"
sample=rbind(TCGA_smaple,GEO_sample)

#ComBat去除批次效应
mod = model.matrix(~1,data = sample)
batch=sample$batch
data0_combat=ComBat(dat=data0,batch = batch,mod=mod,par.prior = T)

#绘制去批次前箱线图及PCA图
##箱线图
pdf(file = "box_FPKM_to_GEO.pdf")
boxplot(log2(data0+1),outline=T, notch=T,las=2)
dev.off()
##PCA图
batch=sample$batch
dat0=as.data.frame(t(data0))
dat.pca0 <- PCA(dat0, graph = FALSE)
pca_plot0 <- fviz_pca_ind(dat.pca0,geom.ind = "point",
                          col.ind = batch,
                          palette = c("#00AFBB", "#E7B800"),
                          addEllipses = TRUE, 
                          legend.title = "Groups")
ggsave(plot = pca_plot0,filename ="PCA_FPKM_to_GEO.pdf")

#绘制去批次后箱线图及PCA图
##绘制箱线图
pdf(file = "box_FPKM_to_GEO_ComBat.pdf")
boxplot(log2(data0_combat+1),outline=T, notch=T,las=2)
dev.off()
##绘制PCA图
dat1=as.data.frame(t(data0_combat))
dat.pca0_combat <- PCA(dat1, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca0_combat,geom.ind = "point",
                         col.ind = batch,
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, 
                         legend.title = "Groups")
ggsave(plot = pca_plot,filename ="PCA_FPKM_to_GEO_ComBat.pdf")



