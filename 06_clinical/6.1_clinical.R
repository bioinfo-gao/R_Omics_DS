# 基因表达与临床特征的关联: 卡方检验 + 注释热图 + 分组箱线图
#
# 输入: rocSigExp.txt(来自 04.ROC-Analysis, 前两列为 futime/fustat),
#       clinical.txt(行=样本, 列=临床变量)
# 输出: 每个基因一张 <gene>_heatmap.pdf, 以及每个基因×每个临床变量一张
#       <gene>_clinicalCor_<clinical>.pdf
#
# ⚠ 第 18 行的循环范围有问题, 见该处说明。

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")

#install.packages("ggpubr")

library(limma)
library(ggpubr)
library(ComplexHeatmap)

expFile="rocSigExp.txt"       
cliFile="clinical.txt"       

# setwd("C:/Users/zhen-/Code/R_code/R_For_DS_Omics/06_clinical")

rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

# ⚠ length() 作用于 data.frame 返回列数, 所以 i 会从 1 遍历到最后一列 ——
#   包括第 1、2 列的 futime 和 fustat。它们会被当成「基因」做中位二分,
#   并产出 futime_heatmap.pdf / fustat_heatmap.pdf 这类无意义结果, 且不报错。
#   若输入沿用 04.ROC 的 rocSigExp.txt(前两列是 futime/fustat), 应从 3 开始:
#   for(i in 3:ncol(rt))。此处未自动修改, 因为你的实际输入列序需要确认。
for(i in 1:length(rt[,1:ncol(rt)])){
    
    gene=colnames(rt)[i]
    
    tumorData=as.matrix(rt[gene])
    
    #rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
    
    data=avereps(tumorData)
    
    Type=ifelse(data[,gene]>median(data[,gene]), "High", "Low")
    Type=factor(Type, levels=c("Low","High"))
    
    data1=cbind(as.data.frame(data), Type)
    data1=data1[order(data1[,gene]),] 
    #??ȡ?ٴ??????ļ?
# ⚠ 临床表在【每次循环里都重新读一遍】, 基因多时会重复 IO; 移到循环外更合适
    cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
    cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))
    
    #?ϲ?????
# 按样本名求交集对齐; 注意 data 与 cli 各自取交集后行序一致, 这一步是安全的
    samSample=intersect(row.names(data), row.names(cli))
    data=data[samSample,,drop=F]
    cli=cli[samSample,,drop=F]
    rt1=cbind(data, cli)
    samSample2=intersect(row.names(data1), row.names(cli))
    data1=data1[samSample2,"Type",drop=F]
    cli=cli[samSample2,,drop=F]
    rt2=cbind(data1, cli)
    
    sigVec=c(gene)
    
# 对每个临床变量做卡方检验, 显著性以星号追加到列名上, 供热图注释显示
    for(clinical in colnames(rt2[,2:ncol(rt2)])){
        data=rt2[c("Type", clinical)]
        colnames(data)=c("Type", "clinical")
# 逐变量剔除 "unknow"(原拼写); 注意各变量剔除的样本不同, 故每个检验的 n 可能不同
        data=data[(data[,"clinical"]!="unknow"),]
        tableStat=table(data)
# ⚠ 卡方检验未做多重检验校正; 临床变量一多, 星号里会混入假阳性
        stat=chisq.test(tableStat)
        pvalue=stat$p.value
        Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
        sigVec=c(sigVec, paste0(clinical, Sig))
    }
    
    colnames(rt2)=sigVec
    
    bioCol=c("#0066FF","#FF9900","#FF0000","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
             "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
             "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
             "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
             "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
             "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
    
    colorList=list()
    colorList[[gene]]=c("Low"="blue", "High"="red")
    j=0
    
# ⚠ 这里的循环变量名 cli 与上面的临床数据框 cli 同名, 循环结束后 cli 变成字符串。
#   本脚本因为下一轮开头会重新 read.table 而侥幸无事, 但这是明显的隐患。
    for(cli in colnames(rt2[,2:ncol(rt2)])){
        cliLength=length(levels(factor(rt2[,cli])))
        cliCol=bioCol[(j+1):(j+cliLength)]
        j=j+cliLength
        names(cliCol)=levels(factor(rt2[,cli]))
# ⚠ 无条件给 cliCol 加一个 "unknow" 颜色: 若该变量本来没有 unknow 这一水平,
#   会多出一个取值不存在的图例项。
        cliCol["unknow"]="grey75"
        colorList[[cli]]=cliCol
    }
    
    #??????ͼ
    ha=HeatmapAnnotation(df=rt2, col=colorList)
    zero_row_mat=matrix(nrow=0, ncol=nrow(rt2))
    Hm=Heatmap(zero_row_mat, top_annotation=ha)
    
    #??????ͼ
    pdf(file=paste0(gene,"_heatmap.pdf"), width=7, height=5)
    draw(Hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    dev.off()
    #?ٴ??????Է?????????ͼ?ν???
    for(clinical in colnames(rt1[,2:ncol(rt1)])){
        data=rt1[c(gene, clinical)]
        colnames(data)=c(gene, "clinical")
        data=data[(data[,"clinical"]!="unknow"),]
        #???ñȽ???
        group=levels(factor(data$clinical))
        data$clinical=factor(data$clinical, levels=group)
        comp=combn(group,2)
        my_comparisons=list()
# 注: 这里的 i 与第 18 行外层循环同名, 但 R 的 for 每轮从固定序列取值,
#     外层循环不会被打乱 —— 只是可读性差。
        for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
        #????????ͼ
        boxplot=ggboxplot(data, x="clinical", y=gene, fill="clinical",
                          xlab=clinical,
                          ylab=paste(gene, " expression"),
                          legend.title=clinical)+ 
            stat_compare_means(comparisons = my_comparisons)
        #????ͼƬ
        pdf(file=paste0(gene,"_clinicalCor_", clinical, ".pdf"), width=5.5, height=5)
        print(boxplot)
        dev.off()
    }
}


