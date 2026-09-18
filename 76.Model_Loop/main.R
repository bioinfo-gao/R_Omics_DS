# LASSO 预后模型「随机搜索」循环 —— 反复重跑建模, 只保留达到阈值的那一轮
#
# 输入: data_orign_rt   行=样本, 列 1=futime 2=fustat 3+=基因表达
#       data_orign_cli  行=样本, 列=临床变量 (供独立预后分析)
# 输出: 每成功一轮建一个 "loop_N_result/" 目录, 内含 cvfit/lambda/ROC 三张 pdf
#       与 lasso_geneCoef / lasso_Risk / UniCox / UniSigExp 四个 txt
#
# ⚠ 本函数依赖 5 个全局变量而非参数, 调用前必须先在环境中定义:
#   workspace, mod_AUC_value, loopTime, pFilter, modelpFdlter
#   (modelpFdlter 是原作者拼写, 不是 Filter; 定义成 modelpFilter 会 object not found)

XXD_lasso_MOD_loop <- function(data_orign_rt, data_orign_cli){
setwd(workspace)
mod_AUC=mod_AUC_value
loop=loopTime
  # ⚠ 循环体内 rt/cli 每轮都重置为同一份数据, 没有任何重抽样。
  #   唯一随机性来自下面 cv.glmnet 的折划分(未设 set.seed), 故每轮 lambda.min
  #   与入选基因不同 —— 本质是在反复碰运气找一个好划分。全流程不可复现。
for(z in 1:loop){
  saytimes=paste0("第",z,"次循环")
  tag=paste0("loop_",z)   # 文件/目录名一律用 ASCII, 避免再次产生中文路径
  setwd(workspace)
  rt=data_orign_rt
  cli=data_orign_cli
    UniCox_outTab=data.frame()
    sigGenes=c("futime","fustat")
    # 第一步: 逐基因单因素 Cox 初筛, 保留 p < pFilter 者
    # 3:ncol 按位置取基因列, 依赖 futime/fustat 恰在前两列
    for(i in colnames(rt[,3:ncol(rt)])){
      cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
      coxSummary = summary(cox)
      coxP=coxSummary$coefficients[,"Pr(>|z|)"]
      if(coxP<pFilter){
        sigGenes=c(sigGenes,i)
        UniCox_outTab=rbind(UniCox_outTab,
                            cbind(id=i,
                                  HR=coxSummary$conf.int[,"exp(coef)"],
                                  HR.95L=coxSummary$conf.int[,"lower .95"],
                                  HR.95H=coxSummary$conf.int[,"upper .95"],
                                  pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
        )
      }
    }
    uniSigExp=rt[,sigGenes]
    uniSigExp_out=cbind(id=row.names(uniSigExp),uniSigExp)
    # length() 作用于 data.frame 返回列数而非行数;
    # 故 >=5 的实际含义是「单因素筛完至少还剩 3 个基因」(另 2 列是 futime/fustat)
    if(length(uniSigExp)>=5){
      x=as.matrix(uniSigExp[,c(3:ncol(uniSigExp))])
      y=data.matrix(Surv(uniSigExp$futime,uniSigExp$fustat))
      # ⚠ matrix=1000 不是 glmnet/cv.glmnet 的参数, 疑为 nlambda 或 maxit 之误。
      #   它会被当作 ... 透传, 既不报错也不生效 —— 属静默无效参数。
      fit=glmnet(x, y, family = "cox",matrix=1000)
      cvfit=cv.glmnet(x, y, family="cox",matrix=1000)
      coef=coef(fit, s = cvfit$lambda.min)
      index=which(coef != 0)
      actCoef=coef[index]
      lassoGene=row.names(coef)[index]
      if(length(lassoGene)>=3){
        geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
        trainFinalGeneExp=uniSigExp[,lassoGene]
        myFun=function(x){crossprod(as.numeric(x),actCoef)}
        # riskScore = 入选基因表达 × LASSO 系数 的线性组合(即 lp)
        # lassoGene 与 actCoef 由同一 index 取出, 列序对齐 —— 这里是对的
        trainScore=apply(trainFinalGeneExp,1,myFun)
        outCol=c("futime","fustat",lassoGene)
        # 中位切点取自本批训练数据, 换队列切点就变
        risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
        lasso_outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
        
        # ⚠ times=c(1,3,5) 假定 futime 单位是年; 上游若传进来的是天, 这就是 1/3/5 天
        ROC_lasso=timeROC(T=lasso_outTab$futime,delta=lasso_outTab$fustat,
                          marker=lasso_outTab$riskScore,cause=1,
                          weighting='aalen',
                          times=c(1,3,5),ROC=TRUE)
        # ⚠ 该 AUC 由建模所用的同一批样本算出(in-sample), 天然偏高, 不代表验证性能
        AUC=ROC_lasso[["AUC"]][["t=1"]]
        sameSample=intersect(row.names(cli),row.names(lasso_outTab))
        risk=lasso_outTab[sameSample,]
        cli=cli[sameSample,]
        Tabrisk=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
        uniTab=data.frame()
        for(i in colnames(Tabrisk[,3:ncol(Tabrisk)])){
          cox <- coxph(Surv(futime, fustat) ~ Tabrisk[,i], data = Tabrisk)
          coxSummary = summary(cox)
          uniTab=rbind(uniTab,
                       cbind(id=i,
                             HR=coxSummary$conf.int[,"exp(coef)"],
                             HR.95L=coxSummary$conf.int[,"lower .95"],
                             HR.95H=coxSummary$conf.int[,"upper .95"],
                             pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
          )
        }
        rownames(uniTab)=uniTab$id
        uniTab_risk=uniTab[uniTab[,"id"]=="riskScore",]
        uniTab_risk_p=as.numeric(uniTab_risk$pvalue)
        if(uniTab_risk_p<0.05){
          uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<0.05,]
          rt1=Tabrisk[,c("futime","fustat",as.vector(uniTab[,"id"]))]
          multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
          multiCoxSum=summary(multiCox)
          multiTab=data.frame()
          multiTab=cbind(
            HR=multiCoxSum$conf.int[,"exp(coef)"],
            HR.95L=multiCoxSum$conf.int[,"lower .95"],
            HR.95H=multiCoxSum$conf.int[,"upper .95"],
            pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
          multiTab=as.data.frame(cbind(id=row.names(multiTab),multiTab))
          multiTab_risk=multiTab[multiTab[,"id"]=="riskScore",]
          multiTab_risk_p=as.numeric(multiTab_risk$pvalue)
          # ⚠ 方法学关键: 以下是一个「棘轮」—— 每成功一轮就把 AUC 门槛抬高 0.01,
          #   而循环不断更换随机划分重试, 直到撞上一个 AUC 更高的。
          #   这是在同一批数据上对检验统计量做选择(selection on the test statistic),
          #   产出的 AUC 与 p 值已失去原本含义, 不能作为模型性能的证据。
          #   要得到可信性能, 必须用完全独立的外部队列评估。
          if(AUC<mod_AUC){
            print(paste0(saytimes,"不符合条件，ROC_AUC不符合阈值"),TRUE)}
          if(AUC>mod_AUC&multiTab_risk_p<modelpFdlter){
            mod_AUC=mod_AUC+0.01
            # 注意: 这里 setwd 进结果子目录后本分支内不再切回,
            # 靠下一轮开头的 setwd(workspace) 复位; 中途报错会把后续输出写错位置。
            setwd(workspace)
            dir.create(paste0(tag,"_result"))
            setwd(paste0(tag,"_result"))
            print(paste0(saytimes,"【###符合条件，输出结果，同时递增阈值###】【###符合条件，输出结，同时递增阈值果###】【###符合条件，输出结果，同时递增阈值###】"))
            print(paste0("AUC阈值增加为_",mod_AUC))
            pdf(file=paste0(tag,"_cvfit.pdf"))
            plot(cvfit)
            dev.off()
            pdf(file = paste0(tag,"_lambda.pdf"))
            plot(fit, xvar = "lambda", label = TRUE)
            dev.off()
            write.table(geneCoef,file=paste0(tag,"_lasso_geneCoef.txt"),sep="\t",quote=F,row.names=F)
            write.table(cbind(id=rownames(lasso_outTab),lasso_outTab),file=paste0(tag,"_lasso_Risk.txt"),sep="\t",quote=F,row.names=F)
          
            write.table(UniCox_outTab,file=paste0(tag,"_UniCox.txt"),sep="\t",row.names=F,quote=F)
            write.table(uniSigExp_out,file=paste0(tag,"_UniSigExp.txt"),sep="\t",row.names=F,quote=F)
            unicox_rt=UniCox_outTab[,2:ncol(UniCox_outTab)]
            unicox_rt=as.matrix(unicox_rt)
            rownames(unicox_rt)=UniCox_outTab$id
            dimnames=list(rownames(unicox_rt),colnames(unicox_rt))
            data=matrix(as.numeric(as.matrix(unicox_rt)),nrow=nrow(unicox_rt),dimnames=dimnames)
            data=avereps(data)
            unicox_rt=as.data.frame(data)
            gene <- rownames(unicox_rt)
            hr <- sprintf("%.3f",unicox_rt$"HR")
            hrLow  <- sprintf("%.3f",unicox_rt$"HR.95L")
            hrHigh <- sprintf("%.3f",unicox_rt$"HR.95H")
            Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
            pVal <- ifelse(unicox_rt$pvalue<0.001, "<0.001", sprintf("%.3f", unicox_rt$pvalue))
            bioCol=rainbow(3, s=0.9, v=0.9)
            pdf(file=paste0(tag,"_ROC.pdf"), width=5, height=5)
            plot(ROC_lasso,time=1,col=bioCol[1],title=FALSE,lwd=2)
            plot(ROC_lasso,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
            plot(ROC_lasso,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
            legend('bottomright',
                   c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_lasso$AUC[1])),
                     paste0('AUC at 3 years: ',sprintf("%.03f",ROC_lasso$AUC[2])),
                     paste0('AUC at 5 years: ',sprintf("%.03f",ROC_lasso$AUC[3]))),
                   col=bioCol[1:3], lwd=2, bty = 'n')
            dev.off()
          }
          if(multiTab_risk_p>modelpFdlter){
            print(paste0(saytimes,"不符合条件，多因素独立预后分析风险值无预测价值"))
          }
        }
        if(uniTab_risk_p>0.05){
          print(paste0(saytimes,"不符合条件，单因素独立预后p值不符合阈值"))
        }
        }
        # 与上面第 134 行处是同一判断的重复分支, 逻辑冗余
        if(uniTab_risk_p>0.05){print("不符合条件，单因素独立预后分析风险值无预测价值")}
      }
      if(length(lassoGene)<3){
        print(paste0(saytimes,"不符合条件，lasso回归得到基因过少"))
      }
    }
    if(length(uniSigExp)<5){
      print(paste0(saytimes,"不符合条件，单因素Cox得到基因过少"))
    }
    
  
}

