# 时间依赖 ROC 筛选: 对每个基因算 5 年 AUC, 保留达标者
#
# 输入: surSigExp.txt (来自 03.Survival)
# 输出: ROC.xls, rocSigExp.txt (供 06_clinical / 70.Clinical_Correlation_Screening 使用),
#       以及 <gene>.ROC.pdf
#
# ⚠ 文件末尾那段多时点 ROC 目前是坏的, 见第 65 行; 原作者在第 64 行也已注明。

# ERROR: dependency 'MASS' is not available for package 'doBy'
# ERROR: dependency 'TH.data' is not available for package 'multcomp'
# ERROR: dependencies 'MASS', 'doBy' are not available for package 'pbkrtest'
# installation of package 'TH.data' had non-zero exit status
# installation of package 'multcomp' had non-zero exit status
# installation of package 'future' had non-zero exit status
# installation of package 'Publish' had non-zero exit status
# installation of package 'future.apply' had non-zero exit status

# ==> conda install -c conda-forge r-timeroc -y # great to settle conflict 
# ==> Ctrl + /         comment out a line 
# ==> Ctrl shift K     delete a  line 
# ==> Ctrl + Delete    delete a word 
# ==> Ctrl+K Ctrl+S  Search for the command "Delete All Right" in the Keyboard Shortcuts editor 
# ==> Ctrl+ Alt + x  delete all Right
# ⚠ 下面 5 行 install.packages 未被注释: source 本文件会触发联网安装。缺包时再手动执行。
install.packages("survivalROC")
install.packages("survminer", dependencies = TRUE)
install.packages("doBy", dependencies = TRUE)
install.packages("pbkrtest", dependencies = TRUE)
install.packages("timeROC", dependencies = TRUE)

library(survivalROC)
library(survival)
library(survminer)
library(timeROC)

setwd("C:\\Users\\zhen-\\Code\\R_code\\R_For_DS_Omics\\04.ROC-Analysis")                     

# ⚠ rocFilter=0 意味着 AUC>0 即通过 —— 实际上不做任何过滤, 所有基因都会进 sigGenes
rocFilter=0                                                                  
rt=read.table("surSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)   

outTab=data.frame()

sigGenes=c("futime","fustat")

# ⚠ 这是一句交互式帮助调用, 属调试残留
?timeROC
# 3:ncol 按位置取基因列, 依赖 futime/fustat 恰在前两列
for(i in colnames(rt[,3:ncol(rt)])){
	   roc=timeROC(T=rt$futime, 
	                   delta=rt$fustat, 
	                   marker = rt[,i], 
	                   cause=1, 
	                   weighting='aalen',time=5,ROC=TRUE)
	   if(roc$AUC[2]>rocFilter){
	       sigGenes=c(sigGenes,i)
	       outTab=rbind(outTab,cbind(gene=i,AUC=roc$AUC[2]))
	   }
}

write.table(outTab,file="ROC.xls",sep="\t",row.names=F,quote=F)    

rocSigExp=rt[,sigGenes]

rocSigExp=cbind(id=row.names(rocSigExp),rocSigExp)
rocSigExp
dim(rocSigExp)

write.table(rocSigExp,file="rocSigExp.txt",sep="\t",row.names=F,quote=F)


gene=colnames(rt)[3]
gene
dim(rt)

# this sentence is WRONG, the next script has details 
# ⚠ 这里用的是 rt$time / rt$event, 但 rt 的列名是 futime / fustat,
#   两者都会取到 NULL, timeROC 无法运行 —— 与上一行原注释所说一致。
#   正确写法是 T=rt$futime, delta=rt$fustat。
#   未自动修改: 这段多时点 ROC 的目标基因(colnames(rt)[3])也需要你确认。
ROC_rt=timeROC(T=rt$time, delta=rt$event,
               marker=rt[,gene], cause=1,
               weighting='aalen',
               times=c(1,3,5,10), ROC=TRUE)
ROC_rt

pdf(file=paste0(gene, ".ROC.pdf"), width=5, height=5)
plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=10,col='yellow',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3])),
         paste0('AUC at 10 years: ',sprintf("%.03f",ROC_rt$AUC[4]))),
       col=c("green",'blue','red','yellow'),lwd=2,bty = 'n')
dev.off()
