# 单样本 scRNA-seq 标准流程: 质控 -> 归一化 -> 降维 -> 聚类 -> 细胞注释
#
# 输入: ./filtered_gene_bc_matrices/hg19/ (10X 三件套目录)
# 输出: QC.Rdata; 其余图形均只画到默认设备, 未落盘
#
# ⚠ 三处值得先看: 第 29/44 行(线粒体比例算了却没用来过滤)、
#   第 70 行(单样本跑 Harmony)、第 105 行(写死的 cluster 编号 + 全程无随机种子)。

####GSE239676####
library(stringr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(scRNAtoolVis)
library(clustree)
library(harmony)



# ⚠ rm(list=ls()) 会清空调用者的工作环境; 且它放在 library() 之后,
#   包虽然仍是加载状态, 但此前所有对象都会没。
rm(list = ls())
setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/53.Single_and_Spatial") # setwd("H:/BCAT1/R分析/data/")

#count <- Read10X(data.dir = './GSE239676/', gene.column = 1)  # 提取第1列作为基因名
count <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/") # pbmc

sce <- CreateSeuratObject(
    counts = count,
    project = "GSE239676",  # 这里用具体的字符串
    min.cells = 3,          # 过滤掉在少于3个细胞中表达的基因
    min.features = 200      # 过滤掉表达基因数少于200的细胞（质控）
)



#计算线粒体比例
# ⚠ 这里算了线粒体比例 mt_percent, 但下面第 44 行的过滤条件里【没有它】。
#   线粒体占比是单细胞质控最标准的一条(通常 <10~20%), 漏掉它意味着
#   濒死细胞不会被剔除。小提琴图能看到这一列, 但看到不等于过滤了。
sce@meta.data$mt_percent = PercentageFeatureSet(sce,pattern = "^MT-")
#计算红细胞比例
HB_gene = c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") #定义红细胞基因
HB_m = match(HB_gene,rownames(sce@assays$RNA)) #在seurat中找到红细胞基因索引
HB_genes = rownames(sce@assays$RNA)[HB_m] #得到匹配红细胞基因的行名
HB_genes = HB_genes[!is.na(HB_genes)] #删除NA值（未匹配到的）
#View(sce@meta.data)
# features= 传基因名向量; 若这些基因都不在对象里会报错, 故上面先做了 match + 去 NA
sce@meta.data$HB_percent = PercentageFeatureSet(sce,features = HB_genes)   #有的时候features会报错，可以换成pattern


#单样本
VlnPlot(sce,features = c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent"),
        group.by = "orig.ident",pt.size = 0,ncol = 4)

# ⚠ 注意上一行原注释写的是「这个样本已经处理过, 不可以过滤」, 但紧接着的 subset 仍在过滤 ——
#   注释与代码相互矛盾, 请确认到底该不该过滤。
#单样本过滤  这个样本已经处理过，不可以过滤，否则匹配不上
sce = subset(sce,subset = nFeature_RNA >200 & nFeature_RNA < 5000 &
                 HB_percent < 3 &
                 nCount_RNA < quantile(nCount_RNA,0.97) &
                 nCount_RNA > 1000)

VlnPlot(sce,features = c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent"),
        group.by = "orig.ident",pt.size = 0,ncol = 4)



#table(GSE163558$orig.ident)

save(sce,file = "QC.Rdata")


####标准化归一化以及降维一步到位####
sce = NormalizeData(sce) %>% #数据归一化处理
    FindVariableFeatures(selection.method = "vst",nfeatures = 2500) %>% #寻找高变基因，可以自己修改nfeatures
    ScaleData() %>% #数据标准化
    RunPCA(npcs = 50, verbose = T) #PCA降维，npcs默认为50，设置多运行会更久
DimPlot(sce,reduction = 'pca',group.by = 'orig.ident',raster=FALSE)



####harmony####
library(harmony)
# ⚠ Harmony 是用来消除【多个样本/批次】之间差异的。本文件从头到尾只有一个样本
#   (上面第 39 行原注释也写着「单样本」), group.by.vars='orig.ident' 只有一个水平,
#   此时 Harmony 没有可对齐的批次, 这一步要么无效、要么直接报错。
sce = RunHarmony(sce,group.by.vars = c('orig.ident'))
DimPlot(sce,reduction = 'harmony',group.by = 'orig.ident',raster=FALSE)



## 降维聚类 ##
ElbowPlot(sce,reduction = 'harmony',ndims = 50)
# ⚠ FindNeighbors/FindClusters/RunUMAP 都含随机过程, 而本文件【全程没有 set.seed】。
#   每次重跑得到的 cluster 编号都可能不同 —— 这一点对下面的注释是致命的。
sce = FindNeighbors(sce,reduction = 'harmony',dims = 1:30)
sce = FindClusters(sce,resolution = 0.3)
sce = RunUMAP(sce,reduction = 'harmony',dims = 1:30)

DimPlot(sce,reduction = 'umap',group.by = 'seurat_clusters',label = T,raster = F)+scale_color_d3('category20')
DimPlot(sce,reduction = 'umap',group.by = 'orig.ident',label = F,raster = F)

####
#sce=JoinLayers(sce)
#cell makers
genes_to_check = c(  "MS4A1",'CD79A', #  B cells
                     "CDH5","PECAM1",#Endothelial
                     "EPCAM",'KRT19',#Epithelial
                     "FN1",'DCN','COL1A1',#Fibroblast
                     'CD68',"CD14",#Myeloid cell
                     "RGS5","PDGFA",#Pericyte
                     'MZB1','DERL3',#Plasma cell
                     "CD2","CD3D" # T/NK
)
p = DotPlot(sce,features = genes_to_check,
            assay = "RNA",group.by = "seurat_clusters")+
    theme(axis.text.x = element_text(angle = 90,hjust = 1))+scale_size(range = c(1,6));p


# ⚠⚠ 下面的注释把 cluster 编号【写死】成 5 / 12,14 / 6,7,9,17 ... 
#   这些编号来自某一次特定的聚类结果, 依赖 resolution、Seurat 版本和随机种子。
#   一旦重跑(而上面没有 set.seed), 编号会重排, 于是 B cell 的标签可能贴到内皮细胞上 ——
#   代码照常运行, 不报任何错。
#   稳妥做法: 在 FindClusters 前 set.seed(), 并把注释依据改成 marker 表达而非编号;
#   至少每次重跑后都用上面第 96 行的 DotPlot 重新核对一遍编号再改这张表。
####细胞注释####
a = length(unique(sce$seurat_clusters))-1
celltype = data.frame(ClusterID = 0:a,
                      celltype = "unkown")
celltype[celltype$ClusterID %in% c(5),2] = "B cell"
celltype[celltype$ClusterID %in% c(12,14),2] = "Endothelials"
celltype[celltype$ClusterID %in% c(6,7,9,17),2] = "Epithelials"
celltype[celltype$ClusterID %in% c(13),2] = "Fibroblast"
celltype[celltype$ClusterID %in% c(3,4,10),2] = "Myeloid"
celltype[celltype$ClusterID %in% c(8,16),2] = "Plasma cell"
celltype[celltype$ClusterID %in% c(0,1,2,11),2] = "T/NK cell"
for (i in 1:nrow(celltype)) {
    sce@meta.data[which(sce$seurat_clusters == celltype$ClusterID[i]),"celltype"] = celltype$celltype[i]
}
table(sce$celltype)



# ⚠ 本文件所有图都只是打印到默认设备。若用 source() 跑, 顶层表达式默认不自动打印,
#   图会一张都出不来 —— 需要显式 print(p) 或用 ggsave()。
###ggplot函数作图###
p = DimPlot(sce,reduction = 'umap',group.by = 'celltype',label = T,raster = F)+
    theme(axis.text = element_blank());p

