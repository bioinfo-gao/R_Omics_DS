
library(WGCNA)
library(limma)
library(tidyverse)

rt=read.table("counts.txt",sep="\t",header=T,check.names=F)  
rt=as.matrix(rt)
rownames(rt)=rt[,1] #【将第一列的基因名作为行名】
exp=rt[,2:ncol(rt)] #【提取第二列往后的表达量】
dimnames=list(rownames(exp),colnames(exp)) #【构建一个列表】
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames) #【从上一步的列表中构建一个新的矩阵】
data=avereps(data)
data2=as.data.frame(data) #【将新的矩阵构建为数据框】

#以01A和11A分组，正常放前面，肿瘤放后面
exp_data_T = data2%>% dplyr::select(str_which(colnames(.), "-01A"))#【读取01A的样本（肿瘤）】
nT = ncol(exp_data_T) 
exp_data_N = data2%>% dplyr::select(str_which(colnames(.), "-11A"))#【读取11A的样本（正常）】
nN = ncol(exp_data_N) 
rt= exp_data_T #【仅保留肿瘤组的数据】
colnames(rt)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\", colnames(rt)) #【简化样品名】
data3=as.data.frame(rt)
duplicated(colnames(data3))#【查看是否有重复列】
rt<-data3[,!duplicated(colnames(data3))] #【将重复列删除】
datExpr0 = as.data.frame(t(rt))
gsg = goodSamplesGenes(datExpr0, verbose = 3)#【检查缺缺失值】
gsg$allOK
#【如果有缺失值就在下方删除】
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

###过滤表达量【将表达量过低的基因删除】#####
meanFPKM=0.1  #【设置过滤值】
n=nrow(datExpr0)
datExpr0[n+1,]=apply(datExpr0[c(1:nrow(datExpr0)),],2,mean)#【读取表达量值】
datExpr0=datExpr0[1:n,datExpr0[n+1,] > meanFPKM] #【将表达量值与过滤值进行比较】
dimnames=list(rownames(datExpr0), colnames(datExpr0))#【保存符合条件的值】
data=matrix(as.numeric(as.matrix(datExpr0)), nrow=nrow(datExpr0), dimnames=dimnames)
datExpr0=avereps(data)


###导出过滤后的基因表达矩阵
filtered_fpkm=t(datExpr0)#【将矩阵进行行列转换，转换为基因名在上，样品名在左】
filtered_fpkm=data.frame(rownames(filtered_fpkm),filtered_fpkm)
names(filtered_fpkm)[1]="sample"
head(filtered_fpkm)
#【将矩阵导出】
write.table(filtered_fpkm, file="rt_filter.xls",row.names=F, col.names=T,quote=FALSE,sep="\t")


###样品分群聚类
sampleTree = hclust(dist(datExpr0), method = "average")#【将每个相似的样品进行聚类】
pdf(file = "1_samplenet.pdf", width = 10, height = 10)#【绘制聚类图】
par(cex = 0.5)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)
dev.off()


hlevel=22500#【根据聚类图，查看是否有离群值，并设置删除高度】

###【删除离群值】
pdf(file = "1_samplenet_cut.pdf", width = 10, height = 10)
par(cex = 0.5)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)
abline(h = hlevel, col = "red")#【h】是设置删除高度
dev.off()
###【保留非离群样本】
clust = cutreeStatic(sampleTree, cutHeight = hlevel, minSize = 10)
table(clust)
keepSamples = (clust==1)#【提取非离群数据】
datExpr0 = datExpr0[keepSamples, ] #【构建非离群矩阵】


###载入临床性状数据
traitData = read.table("clinical.txt",row.names=1,header=T,comment.char = "",check.names=F)
dim(traitData)
names(traitData)
fpkmSamples = rownames(datExpr0)
traitSamples =rownames(traitData)
traitRows = match(fpkmSamples, traitSamples)
datTraits = traitData[traitRows,]
rownames(datTraits) 
collectGarbage()

###【绘制非离群样本聚类图】
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
pdf(file="2_sample_map.pdf",width=10,height=10)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",cex.lab = 1,
                    cex.axis = 1, cex.main = 1)
dev.off()


enableWGCNAThreads() #多线程工作
powers = c(1:13) #【无标度拓扑拟合指数】
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file="3_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
###【无标度拓扑拟合指数散点图】
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") #【若没有数值到达0.9则可以降低标准，但最低不低于0.8】
###平均连通性散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

sft #【查看最佳软阈】
softPower =sft$powerEstimate #【最佳软阈】
adjacency = adjacency(datExpr0, power = softPower)


###【进行邻接矩阵转换】
###【转换为TOM矩阵】
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


###【基因构建网络】
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="4_gene.pdf",width=12,height=9)#【绘制图形】
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

###【构建模块识别】
minModuleSize = 25 #【设置模块基因数目，每个模块最大基因数量为5000】
#【构建基因模块】
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)#【查看模块标识】
dynamicColors = labels2colors(dynamicMods) #【设置模块颜色】
table(dynamicColors) #【查看模块颜色】
#【绘制模块图形】
pdf(file="5_Tree.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

###【将相似的模块进行网络构建】
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)#【提取相似模块】
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file="6_module.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()

#【删除离群模块】
MEDissThres = 1.08  #【设置剪切高度】

#【绘制模块的网络图】
pdf(file="6_module_cut.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()

###【绘制模块与性状数据热图】
nGenes = ncol(datExpr0)#【提取模块基因】
nSamples = nrow(datExpr0)#【提取模块样本】
moduleTraitCor = cor(MEs, datTraits, use = "p")#【计算模块相关性】
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)#【计算模块P值】
pdf(file="8_Module.pdf",width=10,height=8)#【绘制性状数据热图】
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
###【导出基因所在的模块】
moduleColors=dynamicColors
probes = colnames(datExpr0)#【读取基因名】
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
geneOrder =order(geneInfo0$moduleColor)#【读取模块】
geneInfo = geneInfo0[geneOrder, ]
#【保存所有模块基因】
write.table(geneInfo, file = "all_genes.txt",sep="\t",row.names=F,quote=F)


###【输出每个模块的基因】
#【使用循环对每个模块的基因进行导出】
for (mod in 1:nrow(table(moduleColors))){  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("module_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}



###【计算MM和GS】
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

###【批量输出性状和模块散点图】
picDir="module_trait"#【创建文件夹】
dir.create(picDir)
for (trait in traitNames){#【对每个性状和模块进行循环绘制图片】
  traitColumn=match(trait,traitNames)  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      pdfFile=paste("9_", trait, "_", module,".pdf",sep="")
      outPdf=paste(picDir,pdfFile,sep="\\")
      pdf(file=outPdf,width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}


###【输出GS_MM数据】
colnames(datExpr0)#【提取基因名】
probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))#【对计算得到的MM和GS数据进行循环】
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))#【对每个模块进行循环】
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.xls",sep="\t",row.names=F)#【导出每个基因对应的模块和MM、GS值】

###【绘制所有表征基因热图】
nGenes = ncol(datExpr0)#【读取基因名】
nSamples = nrow(datExpr0)#【读取样品名】
plotTOM = dissTOM^7
diag(plotTOM) = NA
pdf(file="10_allgenemap.pdf",width=9, height=9)#【绘制热图】
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

###【绘制部分表征基因热图】
nSelect = 400#【设置基因数目】【随机抽取400个基因】
set.seed(10)
select = sample(nGenes, size = nSelect)#【获取样品】
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]#【读取模型】
plotDiss = selectTOM^7
diag(plotDiss) = NA
pdf(file="11_selectgenemap.pdf",width=9, height=9)#【绘制热图】
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()

###【绘制模块特征关系】
pdf(file="12_dendrogram.pdf",width=8, height=8)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()

###【绘制模块特征热图】
pdf(file="13_modulemap.pdf",width=8, height=8)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

###【模块特征基因聚类树热图】
pdf(file="14_dendrogrammap.pdf", width=9, height=14)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
dev.off()



cytoDir="CytoscapeInput"#【构建文件夹】
dir.create(cytoDir)
for (mod in 1:nrow(table(moduleColors))){ #【对每个模块进行循环】
  modules = names(table(moduleColors))[mod]  
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modProbes = probes[inModule]
  modGenes = modProbes  
  modTOM = TOM[inModule, inModule] 
  dimnames(modTOM) = list(modProbes, modProbes)
  edges_File = paste("CytoscapeInput-edges-", modules , ".txt", sep="")
  nodes_File = paste("CytoscapeInput-nodes-", modules, ".txt", sep="")
  outEdge=paste(cytoDir,edges_File,sep="\\")
  outNode=paste(cytoDir,nodes_File,sep="\\")
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = outEdge,
                                 nodeFile = outNode,
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
  }

