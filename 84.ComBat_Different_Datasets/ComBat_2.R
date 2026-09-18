# 合并 TCGA 与 GEO 两个来源并用 ComBat 去批次, 出去批次前后的箱线图与 PCA 对照
#
# 输入: combined_RNAseq_FPKM.txt(TCGA), geo_exp.csv(来自 56 或 63)
# 输出: TCGA/GEO 合并矩阵的四张图 —— box_FPKM_to_GEO(.pdf/_ComBat.pdf)
#       与 PCA_FPKM_to_GEO(.pdf/_ComBat.pdf)
#
# ✅ 本文件比 83.ComBat 做得好的地方: 先各自 normalizeBetweenArrays 再合并,
#   并且画了去批次【前后】的 PCA 对照 —— 这是判断 ComBat 是否过度校正的关键证据。



#引用包
library(limma)
library(sva)
library(FactoMineR)
library(factoextra)
#读取基因表达文件,并对数据进行处理
# ⚠ 与 83 相同: 这一侧未显式设置行名, 依赖 read.table 对「表头错位」文件的自动处理
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
# 只保留两套数据共有的基因; 之后两侧按同一顺序取子集, 故 cbind 的行是对齐的
sameGene=intersect(rownames(data0),rownames(data2))
TCGA=data0[sameGene,,drop=F]
GEO=data2[sameGene,,drop=F]
##TCGA和GEO组间矫正
# 分别做组内分位数标准化, 再合并 —— 避免一方的分布把另一方拉偏
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
# ⚠ 与 83 同样的问题: mod 只含截距。若 TCGA 与 GEO 两批的肿瘤/正常构成比例不同,
#   批次与生物学分组就是部分共线的, 去批次会连真实差异一起削弱。
#   下面的 PCA 前后对照正是用来看这一点: 若去批次后两批完全重叠且组内结构也消失了,
#   多半是校正过头了。
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
# ⚠ PCA() 默认 scale.unit=TRUE(对每个基因标准化), 与 62.GEO_PCA_Scatter_Plot 里
#   prcomp 默认不标准化的做法不同 —— 两张 PCA 图不可直接比较。
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
# 去批次后的 PCA, 与上面第 53 行的那张构成对照
dat1=as.data.frame(t(data0_combat))
dat.pca0_combat <- PCA(dat1, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca0_combat,geom.ind = "point",
                         col.ind = batch,
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, 
                         legend.title = "Groups")
ggsave(plot = pca_plot,filename ="PCA_FPKM_to_GEO_ComBat.pdf")



