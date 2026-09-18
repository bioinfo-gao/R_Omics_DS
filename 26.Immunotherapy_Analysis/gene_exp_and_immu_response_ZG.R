# 设置工作目录
setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/26.免疫治疗分析")

# 加载所需的 R 包
# limma 包用于处理矩阵和去重复样本（avereps函数）
library(limma)
# ggpubr 包用于绘制统计图表和添加统计比较结果（ggviolin, stat_compare_means）
library(ggpubr)

 
# 免疫特征数据文件（通常是免疫检查点、TMB、IPS等） （The Cancer Immunome Atlas, TCIA）
# Based on CTLA-4 和 PD-1 免疫检查点基因的表达水平来对样本进行分类和评分。
# ips_ctla4_neg_pd1_neg低表达/阴性低表达/阴性双阴性 (Double Negative)。肿瘤免疫原性可能较低，对 PD-1 或 CTLA-4 单一疗法反应最差。
# ips_ctla4_neg_pd1_pos低表达/阴性高表达/阳性PD-1 阳性，CTLA-4 阴性。可能对 PD-1/PD-L1 单一抑制剂治疗有较好反应。
# ips_ctla4_pos_pd1_neg高表达/阳性低表达/阴性CTLA-4 阳性，PD-1 阴性。可能对 CTLA-4 单一抑制剂治疗有较好反应。
# ips_ctla4_pos_pd1_pos高表达/阳性高表达/阳性双阳性 (Double Positive)。免疫微环境激活，通常对 PD-1/PD-L1 抑制剂或 联合治疗 反应最好。
tciaFile="TCIA.txt"     
expFile="geneExp.txt"   # 基因表达数据文件

# ----------------- 数据加载与预处理 -----------------

# 读取基因表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
# 假设 rt 的第一列是样本类型（Type），第二列是核心基因的表达量
# 提取第一列的列名作为核心基因名（用于后续图表标题）
gene=colnames(rt)[1]

# 筛选肿瘤样本
# 假设 rt 中有一个名为 "Type" 的列，用于区分 "Tumor" 和其他样本
tumorData=rt[rt$Type=="Tumor", 1, drop=F]
# 转换为矩阵格式
tumorData=as.matrix(tumorData)

# 标准化样本 ID：将 TCGA 样本 ID 格式（如 TCGA-XX-XXXX-01A-XX）简化为前三段（患者ID，如 TCGA-XX-XXXX）
# 这是为了确保多个样本ID对应一个患者时，可以进行平均处理
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))

# 去重复样本：如果一个患者有多个肿瘤样本（通常是初级肿瘤和转移灶），取平均值
data=avereps(tumorData)

# ----------------- 根据中位数对样本进行分组 -----------------

# 根据核心基因的表达量（data[,gene]）中位数进行分组
# Type = "High"：表达量大于中位数
# Type = "Low"：表达量小于或等于中位数
Type=ifelse(data[,gene]>median(data[,gene]), "High", "Low")
# 将分组结果转换为因子，并指定顺序 (Low在前，High在后)
Type=factor(Type, levels=c("Low","High"))
# 将分组结果 Type 添加回数据框 data
data=cbind(as.data.frame(data), Type)

# ----------------- TCIA 数据加载与合并 -----------------

# 读取 TCIA 数据文件（包含免疫相关分数/特征）
ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)

# 合并数据：找到表达数据 (data) 和免疫数据 (ips) 都有的共同样本ID
sameSample=intersect(row.names(ips), row.names(data))
# 根据共同样本ID筛选免疫数据和表达分组数据
ips=ips[sameSample, , drop=F]
data=data[sameSample, "Type", drop=F]
# 将免疫特征数据 (ips) 与基因表达分组 (data$Type) 合并
data=cbind(ips, data)

# ----------------- 设置统计比较参数 -----------------

# 确保 Type 是因子并设置分组顺序
group=levels(factor(data$Type))
data$Type=factor(data$Type, levels=c("Low", "High"))
group=levels(factor(data$Type))

# 准备 ggpubr 所需的比较列表
comp=combn(group, 2) # 计算所有可能的两两组合
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]} # 在本例中，只有 Low vs High 一种组合

# ----------------- 循环绘制小提琴图 -----------------
# ZG 这4张图的P 值 极其显著
# 遍历除最后一列（Type 分组）外的所有列（即所有免疫特征）
for(i in colnames(data)[1:(ncol(data)-1)]){
    # 提取当前免疫特征 (i) 和基因表达分组 (Type)
    rt=data[,c(i, "Type")]
    
    # 绘制小提琴图
    gg1=ggviolin(rt, x="Type", y=i, fill = "Type", # 填充颜色基于分组 Type
                 xlab="", ylab=i,                     # 设置坐标轴标签
                 legend.title=gene,                 # 图例标题设置为核心基因名
                 add = "boxplot", add.params = list(fill="white"))+ # 在小提琴图上叠加白色箱线图
        stat_compare_means(comparisons = my_comparisons) # 添加统计学比较结果（默认Wilcoxon检验）
    
    # 将图表保存为 PDF 文件
    pdf(file=paste0(i, ".pdf"), width=4.8, height=4.25)
    print(gg1)
    dev.off()
}