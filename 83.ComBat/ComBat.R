

#引用包
library(limma)
library(sva)

#读取基因表达文件,并对数据进行处理
rt=read.table("combined_RNAseq_TPM_log.txt", header=T, sep="\t", check.names=F)
data=as.data.frame(rt)
data0=avereps(data)
GTEx=read.table("GTExNormalExp.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(GTEx)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data2=as.data.frame(data)
sameGene=intersect(rownames(data0),rownames(data2))
TCGA=data0[sameGene,,drop=F]
nTCGA = ncol(TCGA) 
GTEx=data2[sameGene,,drop=F]
nGTEx = ncol(GTEx) 
data0=cbind(TCGA,GTEx)
pdf(file = "boxbind.pdf")
boxplot(log2(data0+1),outline=T, notch=T,las=2)
dev.off()
TCGA_smaple=as.data.frame(colnames(TCGA))
colnames(TCGA_smaple)="group"
TCGA_smaple$batch="1"
GTEx_sample=as.data.frame(colnames(GTEx))
colnames(GTEx_sample)="group"
GTEx_sample$batch="2"
sample=rbind(TCGA_smaple,GTEx_sample)
mod = model.matrix(~1,data = sample)
batch=sample$batch
data0_combat=ComBat(dat=data0,batch = batch,mod=mod,par.prior = T)
write.table(data0_combat,file = "TCGA_GTEx_normalize.txt",sep="\t",quote=F)
pdf(file = "boxbind_normal.pdf")
boxplot(log2(data0_combat+1),outline=T, notch=T,las=2)
dev.off()
