
library(limma)            #???冒?
expFile="combined_RNAseq_FPKM.txt"      #?????????募?
geneFile="gene.txt"       #?????斜??募?


#??取?????募??????????萁??写???
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#??取?????斜??募?,??取铜???????鼗????谋???量
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

#????????
outTab=rbind(ID=colnames(geneExp),geneExp)
write.table(outTab, file="铜????????????量.txt", sep="\t", quote=F, col.names=F)

#####??????????
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
RNA=data[,group==0]
conNum=length(group[group==1])       #????????品??目
treatNum=length(group[group==0])     #????????品??目
sampleType=c(rep(1,conNum), rep(2,treatNum))


rt1=rt1=read.table("铜????????????量.txt", header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
cuproptosis=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
cuproptosis=avereps(cuproptosis)
cuproptosis=cuproptosis[rowMeans(cuproptosis)>0,]

#删????????品
group=sapply(strsplit(colnames(cuproptosis),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
cuproptosis=cuproptosis[,group==0]

#?????约???
corFilter=0.3           
pvalueFilter=0.05  
outTab=data.frame()
for(i in row.names(RNA)){
  if(sd(RNA[i,])>0.1){
    test=wilcox.test(data[i,] ~ sampleType)
    if(test$p.value<0.05){
      for(j in row.names(cuproptosis)){
        x=as.numeric(RNA[i,])
        y=as.numeric(cuproptosis[j,])
        corT=cor.test(x,y)
        cor=corT$estimate
        pvalue=corT$p.value
        if((cor>corFilter) & (pvalue<pvalueFilter)){
          outTab=rbind(outTab,cbind(Cuproptosis=j,RNA=i,cor,pvalue,Regulation="postive"))
        }
        if((cor< -corFilter) & (pvalue<pvalueFilter)){
          outTab=rbind(outTab,cbind(Cuproptosis=j,RNA=i,cor,pvalue,Regulation="negative"))
        }
      }
    }
  }
}

#?????????缘慕???
write.table(file="Result.txt",outTab,sep="\t",quote=F,row.names=F)

#??取铜????????RNA?谋???量
cuproptosisRNA=unique(as.vector(outTab[,"RNA"]))
cuproptosisRNAexp=data[cuproptosisRNA,]
cuproptosisRNAexp=rbind(ID=colnames(cuproptosisRNAexp), cuproptosisRNAexp)
write.table(cuproptosisRNAexp,file="cuproptosisExp.txt",sep="\t",quote=F,col.names=F)
