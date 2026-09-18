###加载R包
library(readxl)
library(tidyverse)
library(GEOquery)
library(tidyverse)
library(GEOquery)
library(limma) 
library(affy)
library(stringr)
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
geo_exp_2=as.data.frame(gene_exp_matrix2)
geo_exp_3=as.data.frame(gene_exp_matrix3)
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
group_list3 <- ifelse(str_detect(pdata3$source_name_ch1,"	Breast tumor"), "T",
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

#####进行数据矫正#####
boxplot(bindgeo,outline=T, notch=T,col=talgroup_list, las=2)
dev.off()
bindgeo_normal=normalizeBetweenArrays(bindgeo)
boxplot(bindgeo_normal,outline=T, notch=T,col=talgroup_list, las=2)
range(bindgeo_normal)
bindgeo_normal <- log2(bindgeo_normal+1)
bindgeo_normal <-as.data.frame(bindgeo_normal)
bindgeo_normal=na.omit(bindgeo_normal)
write.csv(bindgeo_normal,file = "bindgeo_exp.csv")
range(bindgeo_normal)
dev.off()

#####进行差异分析#####
design=model.matrix(~talgroup_list)
fit=lmFit(bindgeo_normal,design)#这里要注意，分组的样本与矩阵样本是否相符，不相符则去文件中调整然后读入
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
