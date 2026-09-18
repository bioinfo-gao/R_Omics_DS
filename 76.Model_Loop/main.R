XXD_lasso_MOD_loop <- function(data_orign_rt, data_orign_cli){
setwd(workspace)
mod_AUC=mod_AUC_value
loop=loopTime
for(z in 1:loop){
  saytimes=paste0("第",z,"次循环")
  setwd(workspace)
  rt=data_orign_rt
  cli=data_orign_cli
    UniCox_outTab=data.frame()
    sigGenes=c("futime","fustat")
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
    if(length(uniSigExp)>=5){
      x=as.matrix(uniSigExp[,c(3:ncol(uniSigExp))])
      y=data.matrix(Surv(uniSigExp$futime,uniSigExp$fustat))
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
        trainScore=apply(trainFinalGeneExp,1,myFun)
        outCol=c("futime","fustat",lassoGene)
        risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
        lasso_outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
        
        ROC_lasso=timeROC(T=lasso_outTab$futime,delta=lasso_outTab$fustat,
                          marker=lasso_outTab$riskScore,cause=1,
                          weighting='aalen',
                          times=c(1,3,5),ROC=TRUE)
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
          if(AUC<mod_AUC){
            print(paste0(saytimes,"不符合条件，ROC_AUC不符合阈值"),TRUE)}
          if(AUC>mod_AUC&multiTab_risk_p<modelpFdlter){
            mod_AUC=mod_AUC+0.01
            setwd(workspace)
            dir.create(paste0(saytimes,"_结果"))
            setwd(paste0(saytimes,"_结果"))
            print(paste0(saytimes,"【###符合条件，输出结果，同时递增阈值###】【###符合条件，输出结，同时递增阈值果###】【###符合条件，输出结果，同时递增阈值###】"))
            print(paste0("AUC阈值增加为_",mod_AUC))
            pdf(file=paste0(saytimes,"_cvfit.pdf"))
            plot(cvfit)
            dev.off()
            pdf(file = paste0(saytimes,"_lambda.pdf"))
            plot(fit, xvar = "lambda", label = TRUE)
            dev.off()
            write.table(geneCoef,file=paste0(saytimes,"_lasso_geneCoef.txt"),sep="\t",quote=F,row.names=F)
            write.table(cbind(id=rownames(lasso_outTab),lasso_outTab),file=paste0(saytimes,"_lasso_Risk.txt"),sep="\t",quote=F,row.names=F)
          
            write.table(UniCox_outTab,file=paste0(saytimes,"_UniCox.txt"),sep="\t",row.names=F,quote=F)
            write.table(uniSigExp_out,file=paste0(saytimes,"_UniSigExp.txt"),sep="\t",row.names=F,quote=F)
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
            pdf(file=paste0(saytimes,"_ROC.pdf"), width=5, height=5)
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

