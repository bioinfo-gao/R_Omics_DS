####GSE239676####
library(stringr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(scRNAtoolVis)
library(clustree)
library(harmony)



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
sce@meta.data$mt_percent = PercentageFeatureSet(sce,pattern = "^MT-")
#计算红细胞比例
HB_gene = c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") #定义红细胞基因
HB_m = match(HB_gene,rownames(sce@assays$RNA)) #在seurat中找到红细胞基因索引
HB_genes = rownames(sce@assays$RNA)[HB_m] #得到匹配红细胞基因的行名
HB_genes = HB_genes[!is.na(HB_genes)] #删除NA值（未匹配到的）
#View(sce@meta.data)
sce@meta.data$HB_percent = PercentageFeatureSet(sce,features = HB_genes)   #有的时候features会报错，可以换成pattern


#单样本
VlnPlot(sce,features = c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent"),
        group.by = "orig.ident",pt.size = 0,ncol = 4)

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
sce = RunHarmony(sce,group.by.vars = c('orig.ident'))
DimPlot(sce,reduction = 'harmony',group.by = 'orig.ident',raster=FALSE)



## 降维聚类 ##
ElbowPlot(sce,reduction = 'harmony',ndims = 50)
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



###ggplot函数作图###
p = DimPlot(sce,reduction = 'umap',group.by = 'celltype',label = T,raster = F)+
    theme(axis.text = element_blank());p

