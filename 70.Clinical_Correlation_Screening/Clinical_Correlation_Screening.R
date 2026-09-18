# 基因表达与各临床变量的相关性筛选 (Kruskal-Wallis)
#
# 输入: rocSigExp.txt (前两列为 futime/fustat, 第 3 列起为基因表达),
#       clinical_after.txt (行=样本, 列=临床变量)
# 输出: clinical_Cor.xls —— 每个基因 × 每个临床变量的 p 值, 末列 SigNum 为显著个数
#
# ⚠ 本文件存在一处会静默给出错误 p 值的问题, 见第 29 行 FIXME。


#设置工作目录
# ⚠ setwd("") 会直接报错, 运行前必须填上真实路径
setwd("")  

#读取表达数据文件
exp=read.table("rocSigExp.txt",sep="\t",header=T,check.names=F,row.names=1)     

#读取临床数据文件
cli=read.table("clinical_after.txt",sep="\t",header=T,check.names=F,row.names=1) 

#以65岁为年龄划分
# 年龄二分; 注意 "unknow"(原拼写) 被保留为独立一类, 会作为第三组参与检验
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))

#合并数据
# 丢掉前两列(futime/fustat), 只保留基因表达列
exp=exp[,3:ncol(exp)]
samSample=intersect(row.names(exp),row.names(cli))
exp=exp[samSample,]
cli=cli[samSample,]
pFilter=0.05   

#临床相关性分析，输出表格
outTab=c()
outTab=rbind(outTab,c("id",colnames(cli),"SigNum"))
colnames(outTab)=c("id",colnames(cli),"SigNum")
for(i in colnames(exp)){
  clinicalPvalVector=c()
  sigSum=0
  for(clinical in colnames(cli)){
# FIXME ⚠ cbind 会把数值与字符拼成【字符矩阵】, expression 列因此变成字符串,
#   kruskal.test 随后按【字典序】而非数值大小排秩 —— 不报错、不告警, 但 p 值是错的。
#   实测: expr=c(2,30,400,3,40,100) 分两组, 本写法 p=0.5127, 正确写法 p=0.8273。
#   正确写法: rt1 <- data.frame(expression=exp[,i], clinical=cli[,clinical])
#   本文件未作改动, 由你确认后再改。
    rt1=cbind(expression=exp[,i],clinical=cli[,clinical])
    cliTest<-kruskal.test(expression ~ clinical, data = rt1)
    pValue=cliTest$p.value
    clinicalPvalVector=c(clinicalPvalVector,pValue)
# ⚠ 这里用的是未经多重检验校正的原始 p; 基因数 × 临床变量数 的组合量很大,
#   SigNum 会包含大量假阳性。若要下结论需按 BH 校正后再计数。
    if(pValue<pFilter){
      sigSum=sigSum+1
    }
  }
  geneClinical=c(i,clinicalPvalVector,sigSum)
  outTab=rbind(outTab,geneClinical)
}
write.table(outTab,file="clinical_Cor.xls",sep="\t",col.names=F,row.names=F,quote=F)

