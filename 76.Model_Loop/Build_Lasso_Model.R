# 76.Model_Loop 的驱动脚本: 准备数据 + 设置全局阈值 + 调用 main.R 里的循环函数
#
# 输入: diffExpLevel.txt (行=基因, 列=样本), clinical.txt (行=样本, 含 futime/fustat)
# 输出: 由 main.R 负责, 见该文件头部说明
#
# ⚠ 本文件存在两处会让流程直接跑不起来的问题, 见第 15 行与第 34 行的 FIXME。


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
# FIXME ⚠ 替换串结尾的 "\\" 会给每个列名追加一个字面反斜杠。
#   实测: "TCGA-AB-1234-01A-11R-XXXX" -> "TCGA-AB-1234\" (末尾多一个反斜杠),
#   而去掉结尾的 \\ 才得到期望的 "TCGA-AB-1234"。
#   后果: 下一行 intersect(colnames(rt), rownames(cli)) 匹配到 0 个样本,
#   sameSample 为空 -> rt 变成 0 行 -> 后续全部失效。本文件未作改动, 由你确认后再改。
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
# 单因素 Cox 的 p 阈值(main.R 内使用)
pFilter=0.05 
#设置模型独立预后阈值
# ⚠ 变量名 modelpFdlter 是原作者的拼写(不是 Filter); main.R 里按同样拼写引用,
#   两处必须一致, 改一处会导致 object not found。
modelpFdlter=0.05
#设置模型ROC阈值
mod_AUC_value=0.7
#设置循环次数
loopTime=100
#开始循环
# FIXME ⚠ 仓库里这个文件已不叫 "主代码.R" —— 它长期以乱码名 "╓≈┤·┬δ.R" 存在,
#   现已重命名为 main.R。因此这一行在改名前就已经 source 不到, 现在应改为 source("main.R")。
#   本文件未作改动, 由你确认后再改。
source("主代码.R")
XXD_lasso_MOD_loop(data_orign_rt,data_orign_cli)



