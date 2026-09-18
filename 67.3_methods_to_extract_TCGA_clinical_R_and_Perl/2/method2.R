
#读取所有xml文件
library(XML)
library(tidyverse)
options(stringsAsFactors = F)
TCGA_clinical_xmls <- dir("TCGA_clinical/",
                    pattern = "*.xml$",
                    recursive = T)

#把每一个xml文件转换为矩阵
cldf <- function(x){
  xmlresult <- xmlParse(file.path("TCGA_clinical/",x))
  xmltop2 <- xmlRoot(xmlresult)
  TCGA_clinical <- xmlToDataFrame(xmltop2[2])
  return(t(TCGA_clinical))
}

#临床数据合并
cl <- lapply(TCGA_clinical_xmls,cldf) 
TCGA_cl <- t(do.call(cbind,cl)) 
clinical <- data.frame(TCGA_cl)
#提取样品ID
ID=clinical$bcr_patient_barcode
#提取年龄
age=(clinical$days_to_birth)
#提取性别
gender=clinical$gender
#提取生存时间
time=clinical$days_to_death
#提取生存状态
status=clinical$vital_status
#提取分期
stage_event=clinical$stage_event
#合并信息
TCGA_merge=cbind(ID,
                 age,
                 gender,
                 time,
                 status,
                 stage_event)
#删除缺失值
TCGA_merge[which(TCGA_merge=="")]=NA
TCGA_clinical=na.omit(TCGA_merge)
TCGA_clinical=as.data.frame(TCGA_clinical)
#删除重复ID
duplicated(TCGA_clinical$ID)
TCGA_clinical<-TCGA_clinical[!duplicated(TCGA_clinical$ID),]
#导出文件
rownames(TCGA_clinical)=TCGA_clinical$ID
TCGA_clinical=TCGA_clinical[,2:ncol(TCGA_clinical)]
write.csv(TCGA_clinical,file = "TCGA_merge.csv",quote = F)
TCGA_clinical=read.csv("TCGA_merge.csv")
TCGA_clinical$age=abs(TCGA_clinical$age)/365
write.csv(TCGA_clinical,file = "TCGA_clinical.csv",quote = F)
