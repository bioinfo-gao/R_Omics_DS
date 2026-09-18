# 环境搭建速记: 安装本目录分析所需的 CRAN / Bioconductor 包
#
# ⚠ 这是一份【随手记】而不是可直接 source 的脚本, 至少有四处跑不通:
#   第 3 行  install.packages("scran") —— scran 在 Bioconductor, 不在 CRAN;
#   第 26/30 行 foobarbaz 是官方文档里的【占位包名】, 并不存在;
#   第 32 行 library() 一次只接受一个包, 第 2、3 个参数其实是 help 与 pos;
#   第 33 行 getSRAfile 是 SRAdb 里的【函数】, 不是包, library() 会报错。
#   另外第 36 行的 sra_con 在本文件中从未创建。
# 建议按需逐行手动执行, 不要整文件 source。

install.packages("Rtsne")
install.packages("RCurl") # library("RCurl") #sudo apt-get install libcurl4-gnutls-dev
install.packages("scran")
install.packages("pheatmap")
install.packages("UpSetR")
install.packages("VennDiagram")

library(UpSetR)

if (!require("BiocManager", quietly = TRUE))  install.packages("BiocManager")

# There's a trick to this where one needs to add biocViews: to the package Description. 
# That's the only solution I've ever seen to allowing automatic installation of bioconductor dependencies.
BiocManager::install("GenomicRanges")
BiocManager::install("SummarizedExperiment")
BiocManager::install("DESeq2")
BiocManager::install("SRAdb")

library(SRAdb)

help(SRAdb) # BiocManager::install("getSRAfile") #https://www.biostars.org/p/93494/


# https://www.biostars.org/p/93494/
library(remotes)
install_version("foobarbaz", "0.1.2")
#An alternative is to install from the GitHub CRAN mirror.

library(remotes)
install_github("cran/foobarbaz")

library("GenomicRanges","SummarizedExperiment","DESeq2")
library("getSRAfile")

# SRP133642 -O /home/gao/Desktop/Code/Bioinfo/scRNA/
getSRAfile(in_acc = c("SRP133642"), sra_con = sra_con,
           destDir = "/home/gao/Desktop/Code/Bioinfo/scRNA/", fileType = 'sra', srcType='ftp')