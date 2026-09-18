# edgeR 差异分析 (TCGA 肿瘤 vs 癌旁) —— TCGA_clinical 系列的第 5 步
#
# 输入: sampleExp.txt (行=基因, 第 1 列基因名, 其余为样本)
# 输出: edgerOut.xls / diffSig.xls / up.xls / down.xls /
#       normalizeExp.txt / diffmRNAExp.txt / vol.pdf
#
# 注: 本文件来自第三方教学材料(文件头尾的 biowolf 链接), 非本项目自写。
# ⚠ 第 23 行的分组是写死的样本数, 见该处说明。

###Video source: http://study.163.com/u/biowolf
######速科生物: http://www.biowolf.cn/
######作者QQ：2749657388

#source("http://bioconductor.org/biocLite.R")   #source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")

# ⚠ 变量名叫 foldChange, 但下面是拿它与 logFC 比较, 所以 2 的含义是 4 倍而非 2 倍
foldChange=2
padj=0.05

# ⚠ setwd 指向的是教程作者的机器路径(lexb4), 在你这里必然不存在, 运行前需改
setwd("C:\\Users\\lexb4\\Desktop\\TCGAclinical\\05.diff")                    #设置工作目录
library("edgeR")
rt=read.table("sampleExp.txt",sep="\t",header=T,check.names=F)  #改成自己的文件名
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
# counts 上的弱过滤
data=data[rowMeans(data)>1,]

#group=c("normal","tumor","tumor","normal","tumor")
# ⚠ 分组是【按样本数写死的】: 前 49 个当正常、后 123 个当肿瘤, 共 172 个。
#   它与 data 的列名没有任何对应关系, 完全依赖 sampleExp.txt 的列恰好按此顺序排列。
#   换一批数据、少一个样本、或列序变化, 标签就整体错位 —— 长度仍相等时不会报错。
#   稳妥做法: 从 TCGA barcode 第 4 段派生分组(本仓库多数脚本都是这么做的)。
group=c(rep("normal",49),rep("tumor",123))                         #按照癌症和正常样品数目修改
# ⚠ design 在本文件中从未被使用 —— 下面走的是 exactTest 路线, 不需要设计矩阵
design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("normal","tumor"))
topTags(et)
ordered_tags <- topTags(et, n=100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
# pseudo.counts 是 edgeR 在 common dispersion 下反算的等效计数, 只适合可视化
newData=y$pseudo.counts

write.table(diff,file="edgerOut.xls",sep="\t",quote=F)
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)   #输出所有基因校正后的表达值（normalizeExp.txt）
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)         #输出差异基因校正后的表达值（diffmRNAExp.txt）

#volcano
# ⚠ 这张火山图的两个轴是【反的】: 常规火山图 x 轴放 logFC、y 轴放 -log10(p),
#   这里恰好对调, 图形相当于旋转了 90 度。不算错, 但与文献惯例不同, 读图时留意。
# ⚠ 下面 yMax=12 是写死的: |logFC| 超过 12 的基因会被画到坐标范围之外而看不见。
pdf(file="vol.pdf")
xMax=max(-log10(allDiff$FDR))+1
#xMax=100
yMax=12
plot(-log10(allDiff$FDR), allDiff$logFC, xlab="-log10(FDR)",ylab="logFC",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC>foldChange,]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub=allDiff[allDiff$FDR<padj & allDiff$logFC<(-foldChange),]
points(-log10(diffSub$FDR), diffSub$logFC, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()

###Video source: http://study.163.com/u/biowolf
######速科生物: http://www.biowolf.cn/
######作者QQ：2749657388

