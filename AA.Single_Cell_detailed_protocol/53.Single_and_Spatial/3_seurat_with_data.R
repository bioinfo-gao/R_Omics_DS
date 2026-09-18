# Seurat 参考映射示例: 用 panc8 数据集做「整合 -> 标签迁移 -> 投影」全流程
#
# 输入: SeuratData 的 panc8 数据集(需先 InstallData("panc8"))
# 输出: pan.UMAO.pdf / violin1.pdf / pan.UMAP3.pdf
#
# 注: 这是 Seurat 官方 vignette 的改写, 用于学习参考映射流程, 非本项目数据。
# ⚠ 第 29 行起的三处 pdf() 里都缺 print(), 见该处说明。

#BiocManager::install("SingleCellExperiment")
#BiocManager::install("scater")
# BiocManager::install("seurat")
# install.packages('Seurat')
#devtools::install_github('satijalab/seurat-data')
library(Seurat)
library(SeuratData)
library(ggplot2)

#InstallData("panc8")
panc8 <- LoadData("panc8")
table(panc8$tech)

# 用两种技术平台的细胞作参考集, 另两种作查询集(第 45 行), 模拟跨平台标签迁移
pancreas.ref <- subset(panc8, tech %in% c("celseq2", "smartseq2"))
# Seurat v5 写法: 按 tech 把 RNA assay 拆成多个 layer, 供后面 IntegrateLayers 使用
pancreas.ref[["RNA"]] <- split(pancreas.ref[["RNA"]], f = pancreas.ref$tech)

# pre-process dataset (without integration)
pancreas.ref <- NormalizeData(pancreas.ref)
pancreas.ref <- FindVariableFeatures(pancreas.ref)
pancreas.ref <- ScaleData(pancreas.ref)
pancreas.ref <- RunPCA(pancreas.ref)
pancreas.ref <- FindNeighbors(pancreas.ref, dims = 1:30)
pancreas.ref <- FindClusters(pancreas.ref)

pancreas.ref <- RunUMAP(pancreas.ref, dims = 1:30)

setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/53.Single_and_Spatial")

# ⚠ 这三处 pdf() ... dev.off() 之间的绘图表达式都没有用 print() 包住。
#   交互式运行时靠自动打印能出图, 但用 source() 执行时顶层表达式默认不打印,
#   结果是生成三个【空白 PDF】且不报错。改成 print(DimPlot(...)) 即可。
pdf(file=paste0("pan.UMAO.pdf"), width=5, height=4.5)

DimPlot(pancreas.ref, group.by = c("celltype", "tech"))   # 已修正拼写: 原为 celltytpe

dev.off()

## ================================================================
## ================================================================


# 至此才做整合; 上面第 21-25 行是未整合状态的对照, 两者可比较批次效应
pancreas.ref <- IntegrateLayers(object = pancreas.ref, method = CCAIntegration, orig.reduction = "pca",
                                new.reduction = "integrated.cca", verbose = FALSE)
pancreas.ref <- FindNeighbors(pancreas.ref, reduction = "integrated.cca", dims = 1:30)
pancreas.ref <- FindClusters(pancreas.ref)


pancreas.query <- subset(panc8, tech %in% c("fluidigmc1", "celseq"))
pancreas.query <- NormalizeData(pancreas.query)
# ⚠ 锚点是基于未整合的 "pca" 找的, 而下面 MapQuery 用的模型建在 integrated.cca 的 UMAP 上。
#   官方 vignette 就是这么写的, 但两者参照系不同, 解读投影结果时需留意。
pancreas.anchors <- FindTransferAnchors(reference = pancreas.ref, query = pancreas.query, dims = 1:30,
                                        reference.reduction = "pca")

predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.ref$celltype, dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)


# 预测标签与真实标签的一致率 —— 这是本示例唯一的定量评估
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)

table(pancreas.query$predicted.id)
VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")


pdf(file=paste0("violin1.pdf"), width=10, height=6)
VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")
dev.off()

## =================================================
## =================================================

pancreas.ref <- RunUMAP(pancreas.ref, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = pancreas.ref, query = pancreas.query,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(pancreas.ref, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2


pdf(file=paste0("pan.UMAP3.pdf"), width=5, height=4.5)
p1 + p2


dev.off()

