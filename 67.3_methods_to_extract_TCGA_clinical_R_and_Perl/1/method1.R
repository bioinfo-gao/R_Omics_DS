# 方法 1: 用 R 直接读 GDC 下载的 clinical.tsv, 抽取常用临床字段
#
# 输入: ./1/clinical.cart.<日期>/clinical.tsv (GDC Data Portal 导出的临床包)
# 输出: TCGA_clinical.csv
#
# ⚠⚠ 第 36 行的 na.omit 会把所有【仍存活】的患者整行删掉, 见该处说明。
# 注: 同级目录 2/ 与 3_pl_language/ 是另外两种做法(R 与 Perl), 可对照。


library(readr)

setwd("C:/Users/zhen-/Code/R_code/R_Omics_DS/67.3_methods_to_extract_TCGA_clinical_R_and_Perl")

clinical=read_tsv("./1/clinical.cart.2023-02-03/clinical.tsv")

ID=clinical$case_submitter_id

age=clinical$age_at_index

gender=clinical$gender

# ⚠ 只取了 days_to_death。GDC 的规则是: 已故患者填 days_to_death,
#   仍存活的患者该字段为 '-- , 其随访时长在 days_to_last_follow_up 里。
#   做生存分析时两者必须合并成一个 futime, 否则就只剩下事件、没有删失。
time=clinical$days_to_death

status=clinical$vital_status

pathologicT=clinical$ajcc_pathologic_t
pathologicM=clinical$ajcc_pathologic_m
pathologicN=clinical$ajcc_pathologic_n

pathologicStage=clinical$ajcc_pathologic_stage


TCGA_merge=cbind(ID,
                    age,
                    gender,
                    time,
                    status,
                    pathologicT,
                    pathologicM,
                    pathologicN,
                    pathologicStage)

# '-- 是 GDC 表示缺失的占位符, 这里统一转成 NA
TCGA_merge[which(TCGA_merge=="'--")]=NA
# ⚠⚠ 承上: 由于 time 列对所有存活患者都是 NA, 这一行 na.omit 会把他们【整批删除】,
#   最终 TCGA_clinical.csv 里只剩已故患者。
#   用这样的队列做 KM / Cox, 相当于没有删失数据, 生存率会被严重低估, 且不报任何错。
#   正确做法: time = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death),
#   再配合 status 区分事件与删失; 之后只对确实必需的列做 na.omit。
#   未自动修改: 需要你确认 clinical.tsv 里随访字段的确切列名。
TCGA_clinical=na.omit(TCGA_merge)
TCGA_clinical=as.data.frame(TCGA_clinical)

# ⚠ 这一行只是把结果打印出来, 真正去重的是下一行
duplicated(TCGA_clinical$ID)
TCGA_clinical<-TCGA_clinical[!duplicated(TCGA_clinical$ID),]


rownames(TCGA_clinical)=TCGA_clinical$ID
TCGA_clinical=TCGA_clinical[,2:ncol(TCGA_clinical)]
write.csv(TCGA_clinical,file = "TCGA_clinical.csv",quote = F)
