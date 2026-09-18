
library(survival)
library("glmnet")
library(survival)
library(survminer)
library(timeROC)
library(limma)
library(future.apply)
library(parallel)
library(caret)
#读取数据
rt=read.table("diffExpLevel.txt",header=T,sep="\t",check.names=F,row.names=1)  
cli=read.table("clinical.txt", header=T, sep="\t", check.names=F, row.names=1)
cli$futime=cli$futime/365
colnames(rt)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\", colnames(rt))
sameSample=intersect(colnames(rt),rownames(cli))
rt=t(rt)
rt=rt[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(cli[,1:2],rt)
data_orign_rt=rt
data_orign_cli=cli
workspace=getwd()
setwd(workspace)
#设置阈值
pFilter=0.05 
#设置模型独立预后阈值
modelpFdlter=0.05
#设置模型ROC阈值
mod_AUC_value=0.7
#设置循环次数
loopTime=100
#开始循环
source("主代码.R")
XXD_lasso_MOD_loop(data_orign_rt,data_orign_cli)



