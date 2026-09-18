# 多因素 Cox 比例风险模型 -> risk score -> 中位二分组
# 输入 input.txt 列序契约: 1=样本id 2=futime 3=fustat 4+=基因表达
#   读入后 row.names=1 吃掉第 1 列, 故索引里 1=futime 2=fustat 3+=表达
# 输出 risk.txt, 供下游 KM / timeROC 直接读取

library(survival)

# check.names=F 保住带 "-" 的基因名与 TCGA 样本条码, 不要去掉
rt=read.table("input.txt",header=T,sep="\t",check.names=F,row.names=1)

# 陷阱: 3:ncol 按位置取表达量列, 依赖 futime/fustat 恰好在前两列。
# 列序一变, 生存时间会被静默 log2 掉, 而下一行按列名取值照样不报错。
rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)
rt[,"futime"]=rt[,"futime"]/365          # 假定原始单位为天

# "~ ." 会把所有剩余列当协变量: 任何忘记删掉的临床列/ID 列都会进模型
cox <- coxph(Surv(futime, fustat) ~ ., data = rt)

# 核心概念: step() 按 AIC 双向逐步回归, 不是按 p 值筛。
# 且选模型与算 riskScore 用的是同一批样本 —— 区分度天然乐观,
# 属于 in-sample 拟合优度, 不能当验证结果; 下结论需独立队列或交叉验证。
cox=step(cox,direction = "both")

# type="risk" 返回相对于协变量均值的 hazard ratio, 即 exp(centered lp), 不是线性预测值本身。
# 需要 lp 请用 type="lp"; 跨数据集打分时 centering 基准仍取自训练数据。
riskScore=predict(cox,type="risk",newdata=rt)

# 中位切点由本批样本决定, 换队列切点就变。跨数据集比较必须固定切点而非各自取中位。
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))

# 落盘 id/futime/fustat/riskScore/risk —— 这是本文件对下游的唯一接口
write.table(cbind(id=rownames(cbind(rt[,1:2],riskScore,risk)),cbind(rt[,1:2],riskScore,risk)),file="risk.txt",sep="\t",quote=F,row.names=F)

cox    # 打印入选变量的 coef / HR / p
