
#设置工作目录
setwd("")  

#读取表达数据文件
exp=read.table("rocSigExp.txt",sep="\t",header=T,check.names=F,row.names=1)     

#读取临床数据文件
cli=read.table("clinical_after.txt",sep="\t",header=T,check.names=F,row.names=1) 

#以65岁为年龄划分
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))

#合并数据
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
    rt1=cbind(expression=exp[,i],clinical=cli[,clinical])
    cliTest<-kruskal.test(expression ~ clinical, data = rt1)
    pValue=cliTest$p.value
    clinicalPvalVector=c(clinicalPvalVector,pValue)
    if(pValue<pFilter){
      sigSum=sigSum+1
    }
  }
  geneClinical=c(i,clinicalPvalVector,sigSum)
  outTab=rbind(outTab,geneClinical)
}
write.table(outTab,file="clinical_Cor.xls",sep="\t",col.names=F,row.names=F,quote=F)

