# 三个 GEO 数据集合并后的去批次对比: ComBat vs removeBatchEffect, 再做差异分析
#
# 输入: 联网下载 GSE205185 / GSE29431 / GSE20711 及平台文件(GPL21185, GPL570)
# 输出: 去批次前后的箱线图与 PCA 共 6 张, 以及 deg_all.txt / upanddown.csv
#
# ⚠⚠ 第 150-153 行的分组是【手工写死的样本计数】, 第 206-208 行的「方法选择」
#   实际只有后一行生效。两处都会静默影响结论, 详见各自说明。
# 注: 前半部分与 57.Multi_GEO_Merge_Analysis 高度重复(同样三个数据集、同样的转 ID 流程)。

###加载R包
library(readxl)
library(tidyverse)
library(GEOquery)
library(tidyverse)
library(GEOquery)
library(limma) 
library(affy)
library(stringr)
library(FactoMineR)
library(factoextra)
library(sva)
###下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE205185', destdir=".", AnnotGPL = T, getGPL = T)
class(gset)
gset[[1]]
gset2 = getGEO('GSE29431', destdir=".", AnnotGPL = T, getGPL = T)
class(gset2)
gset2[[1]]
gset3 = getGEO('GSE20711', destdir=".", AnnotGPL = T, getGPL = T)
class(gset3)
gset3[[1]]
#提取子集
plf1<-gset[[1]]@annotation
plf2<-gset2[[1]]@annotation
plf3<-gset3[[1]]@annotation
#提取平台文件
GPL_data<- getGEO(filename ="GPL21185.soft.gz", AnnotGPL = T)
GPL_data_11 <- Table(GPL_data)
# ⚠ GPL_data1 与 GPL_data2 读的是同一个 GPL570 文件, 解析了两遍
GPL_data1<- getGEO(filename ="GPL570.annot.gz", AnnotGPL = T)
GPL_data_22 <- Table(GPL_data1)
GPL_data2<- getGEO(filename ="GPL570.annot.gz", AnnotGPL = T)
GPL_data_33 <- Table(GPL_data2)
#提取表达量
exp <- exprs(gset[[1]])
probe_name<-rownames(exp)
exp2 <- exprs(gset2[[1]])
probe_name2<-rownames(exp2)
exp3 <- exprs(gset3[[1]])
probe_name3<-rownames(exp3)

###############################################
###########                       #############
###########       数据1转ID       #############
###########                       #############
###############################################
# ⚠ match 对表达矩阵中不存在的探针返回 NA, exp[NA,] 产生整行 NA;
#   第 50 行按【基因符号】过滤并不能去掉它们, 要到 na.omit 才清掉。
loc<-match(GPL_data_11[,1],probe_name)
probe_exp<-exp[loc,]
raw_geneid<-(as.matrix(GPL_data_11[,"GENE_SYMBOL"]))
index<-which(!is.na(raw_geneid))
geneid<-raw_geneid[index]
exp_matrix<-probe_exp[index,]
geneidfactor<-factor(geneid)
gene_exp_matrix<-apply(exp_matrix,2,function(x) tapply(x,geneidfactor,mean))
rownames(gene_exp_matrix)<-levels(geneidfactor)
gene_exp_matrix=na.omit(gene_exp_matrix)

###############################################
###########                       #############
###########       数据2转ID       #############
###########                       #############
###############################################
loc2<-match(GPL_data_22[,1],probe_name2)
probe_exp2<-exp2[loc2,]
raw_geneid2<-(as.matrix(GPL_data_22[,"Gene symbol"]))
index2<-which(!is.na(raw_geneid2))
geneid2<-raw_geneid2[index2]
exp_matrix2<-probe_exp2[index2,]
geneidfactor2<-factor(geneid2)
gene_exp_matrix2<-apply(exp_matrix2,2,function(x) tapply(x,geneidfactor2,mean))
rownames(gene_exp_matrix2)<-levels(geneidfactor2)
gene_exp_matrix2=na.omit(gene_exp_matrix2)

###############################################
###########                       #############
###########       数据3转ID       #############
###########                       #############
###############################################
loc3<-match(GPL_data_33[,1],probe_name3)
probe_exp3<-exp3[loc3,]
raw_geneid3<-(as.matrix(GPL_data_33[,"Gene symbol"]))
index3<-which(!is.na(raw_geneid3))
geneid3<-raw_geneid3[index3]
exp_matrix3<-probe_exp3[index3,]
geneidfactor3<-factor(geneid3)
gene_exp_matrix3<-apply(exp_matrix3,2,function(x) tapply(x,geneidfactor3,mean))
rownames(gene_exp_matrix3)<-levels(geneidfactor3)
gene_exp_matrix3=na.omit(gene_exp_matrix3)

#数据结合
geo_exp_1=as.data.frame(gene_exp_matrix)
geo_exp_1=normalizeBetweenArrays(geo_exp_1)#进行组内数据归一化
geo_exp_2=as.data.frame(gene_exp_matrix2)
geo_exp_2=normalizeBetweenArrays(geo_exp_2)#进行组内数据归一化
geo_exp_3=as.data.frame(gene_exp_matrix3)
geo_exp_3=normalizeBetweenArrays(geo_exp_3)#进行组内数据归一化
sameSample=intersect(rownames(geo_exp_1), rownames(geo_exp_2))
sameSample=as.data.frame(sameSample)
sameSample=intersect(rownames(geo_exp_3), sameSample$sameSample)
gene_exp1=geo_exp_1[sameSample,,drop=F]
gene_exp2=geo_exp_2[sameSample,,drop=F]
gene_exp3=geo_exp_3[sameSample,,drop=F]
bindgeo=cbind(gene_exp1,gene_exp2,gene_exp3)

#####读取分组信息#####
##数据1
pdata <- pData(gset[[1]])
group_list <- ifelse(str_detect(pdata$source_name_ch1,"primary breast tumour"), "T",
                     "N")
group_list
group_list = factor(group_list,
                    levels = c("T","N"))
group_list
pdata$group=group_list
##数据2
pdata2 <- pData(gset2[[1]])
group_list2 <- ifelse(str_detect(pdata2$source_name_ch1,"Breast normal tissue from a breast cancer patient"), "N",
                     "T")
group_list2
group_list2 = factor(group_list2,
                    levels = c("N","T"))
group_list2
pdata2$group=group_list2
##数据3
pdata3 <- pData(gset3[[1]])
group_list3 <- ifelse(str_detect(pdata3$source_name_ch1,"Breast tumor"), "T",
                      "N")
group_list3
group_list3 = factor(group_list3,
                     levels = c("T","N"))
group_list3
pdata3$group=group_list3

#####分组信息合并#####
group1<-(as.matrix(pdata[,"group"]))
row.names(group1)=rownames(pdata)
colnames(group1)="group"
group2<-(as.matrix(pdata2[,"group"]))
row.names(group2)=rownames(pdata2)
colnames(group2)="group"
group3<-(as.matrix(pdata3[,"group"]))
row.names(group3)=rownames(pdata3)
colnames(group3)="group"
talgroup=as.data.frame(rbind(group1,group2,group3))
talgroup_list=factor(talgroup$group,levels = c("N","T"))
write.csv(talgroup,file = "group.csv")

#####多个数据去批次#####
##处理分组
# ⚠⚠ 这两组向量是【按样本数手工写死】的:
#   batchType 假定三批依次是 22 / 66 / 90 个样本(共 178);
#   modType 则是手抄的 T/N 序列(12T,3N,2T,1N,3T,1N | 12N,54T | 88T,2N)。
#   它们与 bindgeo 的实际列没有任何名字层面的对应关系, 完全靠顺序与个数吻合。
#   只要 GEO 上游数据更新、任一样本被过滤、或三个数据集的读取顺序改变,
#   标签就会整体错位 —— 长度仍然相等时【不会报错】, ComBat 与 limma 都照跑,
#   得到的是彻底错误的批次与分组。
#   稳妥做法: 像上面的 group_list 那样从各自的 pData() 里派生, 再按列名对齐。
#   未自动修改: 需要你确认三个数据集当前的实际样本数与顺序。
batchType=c(rep(1,22),rep(2,66),rep(3,90))
modType=c(rep("T",12),rep("N",3),rep("T",2),rep("N",1),rep("T",3),rep("N",1),
          rep("N",12),rep("T",54),
          rep("T",88),rep("N",2))
# mod 保留生物学分组(T/N), 使 ComBat 在去批次时不把肿瘤/正常差异一起抹掉 —— 这一点是对的
mod = model.matrix(~modType)
pdf(file = "pre_normal.pdf",width = 10,height = 10)
boxplot(bindgeo,outline=T, notch=T,col=talgroup_list, las=2)#绘制去批次前箱线图
dev.off()
#绘制去批次前PCA图
# 去批次【前】的 PCA, 与第 176 行的那张构成对照
dat.pca <- PCA(as.data.frame(t(bindgeo)), graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
# ⚠ factoextra 文档中 geom.ind 的取值是 "point"(单数), 本文件三处都写成了 "points";
#   同仓库的 84.ComBat_Different_Datasets 用的是单数形式。
#   本机未安装 factoextra, 我无法实测复数写法的实际行为 ——
#   请你跑一次确认这三张 PCA 图里是否真的画出了散点, 若为空图就改成 "point"。
                         geom.ind = "points",#仅显示"点(points)"（但不是“文本(text)”）（show "points" only (but not "text")）
                         col.ind = talgroup_list,
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, 
                         legend.title = "Groups")
pca_plot
ggsave(plot = pca_plot,filename ="prenormal_PCA.pdf")
#####ComBat法#####
# ComBat: 经验贝叶斯估计批次效应, mod 中的生物学变量会被保留
bindgeo_ComBat=ComBat(dat=bindgeo, batch=batchType, #使用ComBat法去批次
                      mod=mod, par.prior=TRUE)

pdf(file = "ComBatnormal.pdf",width = 10,height = 10)
boxplot(bindgeo_ComBat,outline=T, notch=T,col=talgroup_list, las=2)#绘制去批次后箱线图
dev.off()
#绘制去批次后PCA图
dat.pca2 <- PCA(as.data.frame(t(bindgeo_ComBat)), graph = FALSE)
pca_plot2 <- fviz_pca_ind(dat.pca2,
                         geom.ind = "points",#仅显示"点(points)"（但不是“文本(text)”）（show "points" only (but not "text")）
                         col.ind = talgroup_list,
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, 
                         legend.title = "Groups")
pca_plot2
ggsave(plot = pca_plot2,filename ="afternormal_PCA.pdf")

#####removeBatchEffect法#####
# ⚠ removeBatchEffect 未传 design 参数。不传时它只知道批次、不知道要保留什么,
#   会把与批次相关的生物学差异一并回归掉。
#   应写成 removeBatchEffect(bindgeo, batch=batchType, design=mod)。
#   这也是它与上面 ComBat(传了 mod)不可直接对比的原因。
bindgeo_remove=removeBatchEffect(bindgeo,batchType)
pdf(file = "removenormal.pdf",width = 10,height = 10)
boxplot(bindgeo_remove,outline=T, notch=T,col=talgroup_list, las=2)#绘制去批次后箱线图
dev.off()
#绘制去批次后PCA图
dat.pca3 <- PCA(as.data.frame(t(bindgeo_remove)), graph = FALSE)
pca_plot3 <- fviz_pca_ind(dat.pca3,
                          geom.ind = "points",#仅显示"点(points)"（但不是“文本(text)”）（show "points" only (but not "text")）
                          col.ind = talgroup_list,
                          palette = c("#00AFBB", "#E7B800"),
                          addEllipses = TRUE, 
                          legend.title = "Groups")
pca_plot3
ggsave(plot = pca_plot3,filename ="afternormal_PCA_remove.pdf")


# 下面两行是「二选一」的意思, 但在 R 里它们是顺序执行的
#####选择校正后的数据集进行后续的差异分析#####

###按照实际情况进行选择
# ⚠⚠ 第 206 行随即被第 208 行覆盖 —— 无论注释怎么写, 实际生效的永远是
#   最后一次赋值, 即 removeBatchEffect 的结果。若想用 ComBat, 必须把下一行注释掉。
bindgeo_after=bindgeo_ComBat #选择ComBat法

bindgeo_after=bindgeo_remove #选择removeBatchEffect法


#####进行差异分析#####
# ⚠ 差异分析用的是 talgroup_list; 它与 bindgeo_after 的列顺序同样是按位置对应的
design=model.matrix(~talgroup_list)
fit=lmFit(bindgeo_after,design)#这里要注意，分组的样本与矩阵样本是否相符，不相符则去文件中调整然后读入
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
write.table(deg, file = "deg_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
logFC=1
adj.P.Val = 0.05
k1 = (deg$adj.P.Val < adj.P.Val)&(deg$logFC < -logFC)
k2 = (deg$adj.P.Val < adj.P.Val)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
write.csv(deg,file="upanddown.csv")
