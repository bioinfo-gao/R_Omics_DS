# NKT 细胞亚群标志基因 DotPlot 片段
#
# 输入: 当前环境中已存在的 Seurat 对象 NKT(本文件不负责创建, 需先跑上游流程)
# 输出: 绘图到默认设备, 未落盘
#
# 用途: 按功能分组展示 marker —— T 样 / NK 样 / 活化 / 耗竭四类。
# ⚠ features 用 list 传入时, DotPlot 会按分组加分隔条; 前提是这些基因都在
#   NKT 的默认 assay 里, 缺失的基因会被静默跳过而不报错, 建议先核对:
#   setdiff(unlist(features), rownames(NKT))

DotPlot(NKT, features = list(
    "T-like NKT" =
        c("CD3E","TRAC", "CD8A", "IL7R","KLF2"),
    "NK-like NKT" =
        c("NKG7", "GNLY", "PRF1","KLRD1","CCL5"),
    "Activated" = c("FOS", "JUN", "CD69"),
    "Exhausted" = c("PDCD1","LAG3","TIGIT")
)) + RotatedAxis()
