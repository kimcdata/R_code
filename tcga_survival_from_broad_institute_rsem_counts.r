###################### TCGA SURVIVAL ANALYSIS FROM BROAD INSTITUTE RSEM NORMALISED RNA-SEQ DATA #######################

loadFirehoseData <- function(datafile, clinicalfile, countthreshold = 16, sample.code = "01", pheno.filter = FALSE, pheno){
  
  # pheno is used to search for variable / value pairs within the clinical data. Example: c("metastatic stage","M1")
  
  
  read.delim(clinicalfile,row.names = 1,F,stringsAsFactors = FALSE) -> clin
  t(clin) -> clin
  rownames(clin) <- toupper(clin[,c("patient.bcr_patient_barcode")])
  #

  read.table(datafile,sep="\t",row.names = 1,header = TRUE,stringsAsFactors = FALSE) -> dat.raw
  dat.raw[-1,] -> dat
  
  dimnames(dat) -> dat.dim
  as.matrix(dat) -> dat
  apply(dat, 2, as.numeric) -> dat
  dimnames(dat) <- dat.dim
  
  sample.types <- sapply(strsplit(colnames(dat),"\\."), function(x)x[4])
  
  print("sample types:")
  print(table(sample.types))
  
  colnames(dat) <- sapply(strsplit(colnames(dat),"\\."), function(x)paste(x[1],x[2],x[3],sep="-"))
  rownames(dat) <- sapply(strsplit(rownames(dat),"\\|"),function(x)x[1])
  dat[-c(1:29),] -> dat
  
  dat.tum <- dat[,grep("01",sample.types)]
  
  intersect(colnames(dat.tum),rownames(clin)) -> master.sample.list
  dat.tum[,master.sample.list] -> dat.tum
  clin[master.sample.list,] -> clin.tum
  apply(dat.tum, 1, mean) -> dat.mean
  which(dat.mean > 16) -> sel
  dat.tum[sel,] -> dat.sel
  dat.sel[which(dat.sel < 1)] <- 1
  log2(dat.sel) -> tum.log
  
  by(tum.log, rownames(tum.log), function(x)colMeans(x)) -> tum.by
  t(sapply(tum.by, function(x)x)) -> tum.u
  apply(tum.u, 2, unlist) -> tum.log
  
  
  #print(length(pheno))
  if(pheno.filter){
    if(length(pheno) > 2){
      if(length(pheno) %% 2 ==0){
        print("multi")
        terms <- matrix(pheno, ncol=2, byrow = T)
        hits <- apply(terms, 1, function(x)grep(x[2],clin.tum[,x[1]]))
        Reduce(intersect, hits) -> hits
        print("hits")
        print(hits)
        tum.log[,hits] -> tum.log
        clin.tum[hits,] -> clin.tum
      }
      else {
        print("pheno not passed in pairs")
      }
    } 
    else {
      print("pair")
      hits <- grep(pheno[2],clin.tum[,pheno[1]])
      print("hits:")
      print(hits)
      
      tum.log[,hits] -> tum.log
      clin.tum[hits,] -> clin.tum
    }
  }
  
  tum.log = data.frame(tum.log, stringsAsFactors = F)
  clin.tum = data.frame(clin.tum, stringsAsFactors = F)
  
  return(list(expr = tum.log, clin = clin.tum))
}





survival_with_covar <- function(stime, event, expr, covariates=list(), split, smin=10, gene.label="", data.label="Gene Expression", boxplot.ylab="Gene Expression", include.covar=F,plotKM=FALSE){
  
  # if(!is.list(covariates)){
  # return("Error - X isn't a list")
  # }
  
  require(survival)
  
  if(length(covariates)>0 && include.covar){
    attach(covariates)
  }

  n_ind <- length(expr)
  
  #event <- sapply(event, function(x)if(tolower(x)=="dead"){return(1)}else{return(0)})
  surv.obj <- Surv(as.numeric(stime), as.numeric(event))
  
  beta <- coxph(surv.obj ~ expr)$coefficients[1]
  beta.order <- order(beta*expr, decreasing=T)
  
  sapply(1:n_ind, function(x){
    if(x < smin){
      return(1)
    }
    else if(x > n_ind-smin){
      return(1)
    }
    else {
      fac <- rep(0, n_ind)
      fac[beta.order[1:x]] <- 0
      fac[beta.order[x+1:n_ind]] <- 1
      #summary(coxph(rnd3.surv~factor(rnd3.fac)+as.numeric(rnd3.tab$age.at.diagnosis)))$logtest[3] -> p
      summary(coxph(surv.obj~factor(fac), method="breslow"))$coefficients[,"Pr(>|z|)"][1] -> p
      return(p)
    }
  }) -> p.opti
  
  which.min(p.opti) -> cutoff
  
  fac <- rep(0, n_ind)
  fac[beta.order[1:cutoff]] <- "bHigh Risk"
  fac[beta.order[cutoff+1:n_ind]] <- "aLow Risk"
  fac <- as.factor(fac)
  
  if(include.covar){
    on.exit(detach(covariates))
    form <- as.formula(paste("surv.obj ~ fac+",paste(names(covariates),collapse="+")))
  } else {
    form <- as.formula("surv.obj ~ fac")
  }
  
  coxph(form, method="breslow") -> cox.model
  #coxph(surv.obj ~ factor(fac), method="breslow") -> cox.model
  #coxph(surv.obj ~ factor(fac)*age*gender, method="breslow") -> cox.model.int
  #eval(parse(text=paste("cox.model.int <<- coxph(surv.obj~",paste(colnames(covar),collapse="*"),",method=\"breslow\",data=covar)")))
  survfit(surv.obj ~ factor(fac)) -> sf.model
  #survdiff(surv.obj ~ factor(fac)) -> sd.model
  
  summary(cox.model) -> cox.summ
  #print(cox.summ)
  
  round(cox.summ$coefficients[1,"Pr(>|z|)"],6) -> pv
  round(cox.summ$conf.int[1,"exp(coef)"],2) -> hr
  round(cox.summ$conf.int[1,"lower .95"],2) -> ci.l
  round(cox.summ$conf.int[1,"upper .95"],2) -> ci.u
  
  #x11();
  #par(mfrow=c(1,2))
  if(plotKM){
    pdf(file=paste(gene.label,".",data.label,".pdf",sep=""), width=12, height=7)
    nf <- layout(matrix(c(1,2), ncol=2), widths=c(4,2))
    #layout.show(nf)
    plot(sf.model, col=c("blue", "red"), xlab="Survival - Days", ylab="Probability of survival", main=paste("Survival analysis - ",data.label,": ",gene.label,". \n HR = ",hr,". CI(L) = ",ci.l,". CI(U) = ",ci.u,". P-value = ",pv,sep=""), las=1)
    legend("topright", col=c("red", "blue"), legend=c(paste("High Risk: ",cutoff,sep=""), paste("Low Risk: ", n_ind-cutoff,sep="")), pch="---")
    print("High risk group")
    print(cutoff)
    print("Low risk group")
    print(n_ind-cutoff)
    
    expr.p <- t.test(expr[fac=="bHigh Risk"], expr[fac=="aLow Risk"])$p.value
    #x11(); 
    boxplot(expr~fac, main=paste("Expression: ", gene.label, "\n p-value = ",round(expr.p, 6), sep=""), las=1, ylab=boxplot.ylab, xaxt="none") -> bx
    axis(1, at=c(1,2), labels=c("Low Risk","High Risk"))
    
    graphics.off()
  }
  
  return(list(pv=pv, hr=hr))
  
}


collapseSurvData <- function(ttd, ttlf, vitalStatus){
  
  survData <- data.frame(ttlf = ttlf, ttd = ttd, vs = vitalStatus, stringsAsFactors = F)
  
  apply(survData, 1, function(x){
    
    if(x[3] == "alive"){
      return(c(stime = x[1], event = 0))
    } else {
      return(c(stime = x[2], event = 1))
    }
    
  }) -> survData
  
  t(survData) -> survData
  
  survTimeIsInt = !is.na(survData[,1])
  
  survData = data.frame(survData, isNA = survTimeIsInt, stringsAsFactors = F) 
  
  return(survData)
  
}




