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
batchType=c(rep(1,22),rep(2,66),rep(3,90))
modType=c(rep("T",12),rep("N",3),rep("T",2),rep("N",1),rep("T",3),rep("N",1),
          rep("N",12),rep("T",54),
          rep("T",88),rep("N",2))
mod = model.matrix(~modType)
pdf(file = "pre_normal.pdf",width = 10,height = 10)
boxplot(bindgeo,outline=T, notch=T,col=talgroup_list, las=2)#绘制去批次前箱线图
dev.off()
#绘制去批次前PCA图
dat.pca <- PCA(as.data.frame(t(bindgeo)), graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "points",#仅显示"点(points)"（但不是“文本(text)”）（show "points" only (but not "text")）
                         col.ind = talgroup_list,
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, 
                         legend.title = "Groups")
pca_plot
ggsave(plot = pca_plot,filename ="prenormal_PCA.pdf")
#####ComBat法#####
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


#####选择校正后的数据集进行后续的差异分析#####

###按照实际情况进行选择
bindgeo_after=bindgeo_ComBat #选择ComBat法

bindgeo_after=bindgeo_remove #选择removeBatchEffect法


#####进行差异分析#####
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
