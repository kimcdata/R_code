######################### SURVIVAL ANALYSIS IN R WITH LARGE GENE EXPRESSION MATRICES #########################

######### INPUTS #########

#survival_times, numeric vector
#event_flag, numeric vector
#gene_expression, numeric vector
#min_group_size, integer
#gene.label, character
#data.type.label, character
#boxplot.ylab, character
#include.covar, boolean
#plotKM, boolean

###################################################################################################################
####################### FETCH SURVFIT OBJECT AFTER FINDING MINIMUM P VALUE ###################################
#######################################################################################################################


gex_survival_getfit <- function(survival_time, event_flag, gene_expression, covariates=list(), min_group_size=10, include.covar=F){
  
  require(survival)
  
  if(length(covariates)>0 && include.covar){
    attach(covariates)
  }
  if(!is.numeric(survival_time) || !is.numeric(event_flag) || !is.numeric(gene_expression)){
    print("survival_time, event & gene_expression must be numeric")
    return(0)
  }
  
  n_ind = length(gene_expression)
  
  surv.obj = Surv(survival_time, event_flag)
  
  beta = coxph(surv.obj ~ gene_expression)$coefficients[1]
  beta.order = order(beta*gene_expression, decreasing=T)
  
  sapply(1:n_ind, function(x){
    if(x < min_group_size){
      return(1)
    }
    else if(x > n_ind-min_group_size){
      return(1)
    }
    else {
      fac = rep(0, n_ind)
      fac[beta.order[1:x]] = 0
      fac[beta.order[x+1:n_ind]] = 1
      p = summary(coxph(surv.obj~factor(fac), method="breslow"))$coefficients[,"Pr(>|z|)"][1]
      return(p)
    }
  }) -> p.opti
  
  which.min(p.opti) -> cutoff
  
  fac = rep(0, n_ind)
  fac[beta.order[1:cutoff]] = "bHigh Risk"
  fac[beta.order[cutoff+1:n_ind]] = "aLow Risk"
  fac = as.factor(fac)
  
  if(include.covar){
    on.exit(detach(covariates))
    form = as.formula(paste("surv.obj ~ fac+",paste(names(covariates),collapse="+")))
  } else {
    form = as.formula("surv.obj ~ fac")
  }
  
  sf.model = survfit(form)
  
  return(sf.model)
  
}

####################################################################################################################
############################# FETCH PATIENT SPLIT WITH SMALLEST P VALUE ############################################
####################################################################################################################



gex_survival_getsplit <- function(survival_time, event_flag, gene_expression, min_group_size=10){
  
  require(survival)
  
  if(!is.numeric(survival_time) || !is.numeric(event_flag) || !is.numeric(gene_expression)){
    print("survival_time, event & gene_expression must be numeric")
    return(0)
  }
  
  n_ind = length(gene_expression)
  
  surv.obj = Surv(survival_time, event_flag)
  
  beta = coxph(surv.obj ~ gene_expression)$coefficients[1]
  beta.order = order(beta*gene_expression, decreasing=T)
  
  sapply(1:n_ind, function(x){
    if(x < min_group_size){
      return(1)
    }
    else if(x > n_ind-min_group_size){
      return(1)
    }
    else {
      fac = rep(0, n_ind)
      fac[beta.order[1:x]] = 0
      fac[beta.order[x+1:n_ind]] = 1
      p = summary(coxph(surv.obj~factor(fac), method="breslow"))$coefficients[,"Pr(>|z|)"][1]
      return(p)
    }
  }) -> p.opti
  
  which.min(p.opti) -> cutoff
  
  fac = rep(0, n_ind)
  fac[beta.order[1:cutoff]] = "High Risk"
  fac[beta.order[cutoff+1:n_ind]] = "Low Risk"
  fac = as.factor(fac)
  
  return(fac)
  
}

###########################################################################################################################################
####################### SURVIVAL ANALYSIS WITH COVARIATES (WITH KM PLOT USING BASE PLOT PACKAGE) ##############################################
##########################################################################################################################################


gex_survival_with_covar <- function(survival_time, event_flag, gene_expression, covariates=list(), min_group_size=10, gene.label="", data.type.label="Gene Expression", boxplot.ylab="Gene Expression", include.covar=F,plotKM=FALSE, file_prefix = ""){
  
  require(survival)
  require(clinfun)
  
  if(length(covariates)>0 && include.covar){
    attach(covariates)
  }
  if(!is.numeric(survival_time) || !is.numeric(event_flag) || !is.numeric(gene_expression)){
    print("survival_time, event & gene_expression must be numeric")
    return(0)
  }
  
  n_ind = length(gene_expression)
  
  surv.obj = Surv(survival_time, event_flag)
  
  beta = coxph(surv.obj ~ gene_expression)$coefficients[1]
  beta.order = order(beta*gene_expression, decreasing=T)
  
  sapply(1:n_ind, function(x){
    if(x < min_group_size){
      return(1)
    }
    else if(x > n_ind-min_group_size){
      return(1)
    }
    else {
      fac = rep(0, n_ind)
      fac[beta.order[1:x]] = 0
      fac[beta.order[x+1:n_ind]] = 1
      p = summary(coxph(surv.obj~factor(fac), method="breslow"))$coefficients[,"Pr(>|z|)"][1]
      return(p)
    }
  }) -> p.opti
  
  which.min(p.opti) -> cutoff
  
  fac = rep(0, n_ind)
  fac[beta.order[1:cutoff]] = "bHigh Risk"
  fac[beta.order[cutoff+1:n_ind]] = "aLow Risk"
  fac = as.factor(fac)
  
  if(include.covar){
    on.exit(detach(covariates))
    form = as.formula(paste("surv.obj ~ fac+",paste(names(covariates),collapse="+")))
  } else {
    form = as.formula("surv.obj ~ fac")
  }
  
  cox.model = coxph(form, method="breslow", x = T)
  sf.model = survfit(form)
  
  cox_err_cv = coxphERR(cox.model, 2:4)
  cox_err_gene = coxphERR(cox.model, 1)
  cox_err_model = coxphERR(cox.model)
  cox_err_diff = cox_err_model["ERR"] - cox_err_gene["ERR"]
  
  
  #print("sf.model")
  #str(sf.model)

  cox.summ = summary(cox.model)
  
  if(include.covar){
    #print(cox.summ)
    pv = cox.summ$coefficients[,"Pr(>|z|)"]
    hr = cox.summ$coefficients[,"exp(coef)"]
    names(pv) = paste("PV-",names(pv),sep="")
    names(hr) = paste("HR-",names(hr),sep="")
  }
  if(!include.covar){
  pv = round(cox.summ$coefficients[1,"Pr(>|z|)"],6)
  hr = round(cox.summ$conf.int[1,"exp(coef)"],2)
  ci.l = round(cox.summ$conf.int[1,"lower .95"],2)
  ci.u = round(cox.summ$conf.int[1,"upper .95"],2)
  }
  
  if(plotKM & !include.covar){
    pdf(file=paste(file_prefix,".",gene.label,".",data.type.label,".pdf",sep=""), width=12, height=7)
    nf <- layout(matrix(c(1,2), ncol=2), widths=c(4,2))

    plot(sf.model, 
         col=c("blue", "red"), 
         xlab="Survival - Days", 
         ylab="Probability of survival", 
         main=paste("Survival analysis - ",
                    data.type.label,": ",
                    gene.label,
                    ". \n HR = ", hr,
                    ". CI(L) = ", ci.l,
                    ". CI(U) = ", ci.u,
                    ". P-value = ",pv,
                    sep=""), 
         las=1)
    
    legend("topright", 
           col=c("red", "blue"), 
           legend=c(paste("High Risk: ",cutoff,sep=""), paste("Low Risk: ", n_ind-cutoff,sep="")), 
           pch="---")
    
    print("High risk group")
    print(cutoff)
    print("Low risk group")
    print(n_ind-cutoff)
    
    bx = boxplot(gene_expression~fac, main=paste("Expression: ", gene.label, sep=""), las=1, ylab=boxplot.ylab, xaxt="none")
    axis(1, at=c(1,2), labels=c("Low Risk","High Risk"))
    
    dev.off()
  }
  
  return(c(cox_err_cv = cox_err_cv, cox_err_model = cox_err_model, cox_err_gene = cox_err_gene, cox_err_diff = cox_err_diff, pv, hr))
  
}

#############################################################################################################################
####################### SURVIVAL ANALYSIS WITH PRESET GROUPS DIVIDING THE DATA ##############################################
#############################################################################################################################


survival_with_covar_preset_groups <- function(survival_times, event_flag, groups, covariates=list(), min_group_size=10, gene.label="", data.type.label="Custom Groups", include.covar=F,plotKM=FALSE){
  
  require(survival)
  
  if(length(covariates)>0 && include.covar){
    attach(covariates)
  }
  if(!is.numeric(survival_times) || !is.numeric(event)){
    print("survival_times, event must be numeric")
    return(0)
  }
  
  n_ind = length(groups)
  
  surv.obj = Surv(survival_time, event_flag)
  
  fac = as.factor(groups)
  
  if(include.covar){
    on.exit(detach(covariates))
    form = as.formula(paste("surv.obj ~ fac+",paste(names(covariates),collapse="+")))
  } else {
    form = as.formula("surv.obj ~ fac")
  }
  
  cox.model = coxph(form, method="breslow")
  sf.model = survfit(surv.obj ~ factor(fac))
  
  cox.summ = summary(cox.model)
  
  pv = round(cox.summ$coefficients[1,"Pr(>|z|)"],6)
  hr = round(cox.summ$conf.int[1,"exp(coef)"],2)
  ci.l = round(cox.summ$conf.int[1,"lower .95"],2)
  ci.u = round(cox.summ$conf.int[1,"upper .95"],2)
  
  if(plotKM){
    pdf(file=paste(data.type.label,".pdf",sep=""), width=12, height=7)

    plot(sf.model, 
         col=c("blue", "red"), 
         xlab="Survival - Days", 
         ylab="Probability of survival", 
         main=paste("Survival analysis - ",
                    data.type.label,": ",
                    ". \n HR = ", hr,
                    ". CI(L) = ", ci.l,
                    ". CI(U) = ", ci.u,
                    ". P-value = ",pv,
                    sep=""), 
         las=1)
    
    legend("topright", 
           col=c("red", "blue"), 
           legend=c(paste("High Risk: ",cutoff,sep=""), paste("Low Risk: ", n_ind-cutoff,sep="")), 
           pch="---")
    
    print("High risk group")
    print(cutoff)
    print("Low risk group")
    print(n_ind-cutoff)
    
    graphics.off()
  }
  
  return(c(pv=pv, hr=hr))
  
}

####################################################################################################################################
####################### CREATING DATAFRAME OF SURVIVAL DATA FROM "Time to death" and "Time to last followup" #######################
####################################################################################################################################



collapse_surv_data <- function(ttd, ttlf, vitalStatus){
  
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



#####################################################################################################################
##################################### coxphERR from library clinfun #################################################
####################################################################################################################

coxphERR_fixed = function (phfit, ngamma = NULL) 
{
  if (class(phfit) != "coxph") 
    stop("phfit shoud be coxph class object")
  if (is.null(phfit$x)) 
    stop("coxph should have been called with x=TRUE option")
  ss <- phfit$n
  covec <- phfit$x - matrix(rep(apply(phfit$x, 2, mean), ss), 
                            nrow = ss, byrow = T)
  if (is.null(ngamma)) 
    ngamma <- 1:ncol(covec)
  covec.gamma <- covec[, ngamma, drop = FALSE]
  covec.beta <- covec[, -ngamma, drop = FALSE]
  xf <- exp(covec %*% phfit$coef)
  xf.gamma <- exp(covec.gamma %*% phfit$coef[ngamma])
  ERR <- (log(mean(xf.gamma)))/(0.577215665 + (log(mean(xf))))
  g <- log(mean(xf))
  h <- log(mean(xf.gamma))
  dgb <- (t(xf) %*% covec.beta)/sum(xf)
  dga <- (t(xf) %*% covec.gamma)/sum(xf)
  dha <- (t(xf.gamma) %*% covec.gamma)/sum(xf.gamma)
  dxfb <- (-h * dgb)/((0.577215665 + g)^2)
  dxfa <- (((0.577215665 + g) * dha) - (h * dga))/((0.577215665 + 
                                                      g)^2)
  dxf <- cbind(dxfa, dxfb)
  newcovmat11 <- phfit$var[ngamma, ngamma, drop = FALSE]
  newcovmat12 <- phfit$var[ngamma, -ngamma, drop = FALSE]
  newcovmat21 <- phfit$var[-ngamma, ngamma, drop = FALSE]
  newcovmat22 <- phfit$var[-ngamma, -ngamma, drop = FALSE]
  newcovmat <- cbind(rbind(newcovmat11, newcovmat21), rbind(newcovmat12, 
                                                            newcovmat22))
  vr <- var(xf.gamma)/ss
  vs <- var(xf)/ss
  cvrs <- cov(xf.gamma, xf)/ss
  dr <- 1/((mean(xf.gamma)) * (0.577215665 + log(mean(xf))))
  ds <- (-log(mean(xf.gamma)))/(mean(xf) * ((0.577215665 + 
                                               log(mean(xf)))^2))
  dd <- cbind(dr, ds)
  newmat <- cbind(rbind(vr, cvrs), rbind(cvrs, vs))
  se.ERR <- sqrt((dd %*% newmat %*% t(dd)) + (dxf %*% newcovmat %*% 
                                                t(dxf)))
  out <- c(ERR, se.ERR)
  names(out) <- c("ERR", "se.ERR")
  out
  
}
