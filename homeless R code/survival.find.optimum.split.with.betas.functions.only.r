getstime <- function(ci){

tab <- cbind(ci$days_to_last_followup, ci$days_to_death, ci$vital_status)

apply(tab, 1, function(x){

if(x[3]=="Dead"){
stime=as.numeric(x[2])
status=1
} else {
stime=as.numeric(x[1])
status=0
}
return(c(stime, status))

}) -> rnd3.surv
return(t(rnd3.surv))

}

survit <- function(stime, event, expr, X=list(), split, smin=(length(stime)/10), gene.label="", data.label="Gene Expression", boxplot.ylab="Gene Expression", include.covar=TRUE,covariates=c(),plotKM=FALSE){

# if(!is.list(covariates)){
# return("Error - X isn't a list")
# }

#attach(X)

n_ind <- length(stime)
surv.obj <- Surv(stime, event)

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

#x11();
#plot(p.opti, type="l", main="Cox model P-values at different splits", ylim=c(0,1), xlab="Split - index", ylab="P value")
#abline(h=0.05, lty="dashed", col="grey20")

#which.min(p.opti) -> pvmin
#pv_col <- rep("black", length(p.opti))
#pv_col[pvmin] <- "red"
#points(1:length(p.opti), p.opti, pch=19, cex=0.5, col=pv_col)

which.min(p.opti) -> cutoff

fac <- rep(0, n_ind)
fac[beta.order[1:cutoff]] <- "bHigh Risk"
fac[beta.order[cutoff+1:n_ind]] <- "aLow Risk"

if(is.data.frame(covariates)){
covar <- data.frame(factor(fac), covariates)
} else {
covar <- data.frame(factor(fac))
}
colnames(covar)[1] <- "x"


coxph(surv.obj ~ factor(fac), method="breslow") -> cox.model
#coxph(surv.obj ~ factor(fac)*age*gender, method="breslow") -> cox.model.int
eval(parse(text=paste("cox.model.int <<- coxph(surv.obj~",paste(colnames(covar),collapse="*"),",method=\"breslow\",data=covar)")))
survfit(surv.obj ~ factor(fac)) -> sf.model
survdiff(surv.obj ~ factor(fac)) -> sd.model

if(include.covar){
summary(cox.model.int) -> cox.summ
} else {
summary(cox.model) -> cox.summ
}
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
plot(sf.model, col=c("blue", "red"), xlab="Survival - Months", ylab="Probability of survival", main=paste("Survival analysis - ",data.label,": ",gene.label,". \n HR = ",hr,". CI(L) = ",ci.l,". CI(U) = ",ci.u,". P-value = ",pv,sep=""), las=1)
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
return(list(C=cox.model, CA = cox.model.int, F=sf.model, D=sd.model))

}

survit(stime=surv.data.sel[,1],event=surv.data.sel[,2],expr=data.rma.sel.unique["RND3",],gene.label="RND3",data.label="REMBRANDT - Expression",boxplot.ylab="Expression",plotKM=TRUE) -> rnd3.survit


survproves <- function(id.sel, t1, t2, smin=(length(id.sel)/10), gene, ageopti=TRUE, plotKM=FALSE, include.covar=TRUE){

id.selc <- intersect(id.sel, rownames(clinical_info.cnv))
id.sele <- intersect(id.sel, rownames(clinical_info.expr))

if(length(which(rownames(cnv_data_combat_unique.sel)==gene))==0){
	print("Gene not found")
	return(list(cnv=NA, expr=NA))
}
if(length(which(rownames(combat_data_unique.sel)==gene))==0){
	print("Gene not found")
	return(list(cnv=NA, expr=NA))
}

surv.obj<- Surv(getstime(na.omit(clinical_info.cnv[id.selc,]))[,1], getstime(na.omit(clinical_info.cnv[id.selc,]))[,2])
age <- clinical_info.cnv[id.selc,"age_at_initial_pathologic_diagnosis"]

if(ageopti){
n_ind <- length(age)
beta <- coxph(surv.obj ~ age)$coefficients[1]
beta.order <- order(beta*age, decreasing=T)

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

fac.age <- rep(0, n_ind)
fac.age[beta.order[1:cutoff]] <- "High Risk - Age"
fac.age[beta.order[cutoff+1:n_ind]] <- "Low Risk - Age"

age <- factor(fac.age)
}

#cova <- data.frame(age=clinical_info.cnv[id.selc,"age_at_initial_pathologic_diagnosis"], gender=clinical_info.cnv[id.selc,"gender"])
cova <- data.frame(age=age)
#x11(width=10);
survit(stime=getstime(na.omit(clinical_info.cnv[id.selc,]))[,1], 
event=getstime(na.omit(clinical_info.cnv[id.selc,]))[,2], 
expr=na.omit(cnv_data_combat_unique.sel[gene, id.selc]), 
X=list(), 
gene.label=gene, 
data.label=t1,
boxplot.ylab="Copy Number Ratio",
smin=smin,
covariates=cova,
plotKM=plotKM,
include.covar=include.covar
) -> cnv


surv.obj<- Surv(getstime(na.omit(clinical_info.expr[id.sele,]))[,1], getstime(na.omit(clinical_info.expr[id.sele,]))[,2])
age <- clinical_info.expr[id.sele,"age_at_initial_pathologic_diagnosis"]

if(ageopti){
n_ind <- length(age)
beta <- coxph(surv.obj ~ age)$coefficients[1]
beta.order <- order(beta*age, decreasing=T)

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

fac.age <- rep(0, n_ind)
fac.age[beta.order[1:cutoff]] <- "High Risk - Age"
fac.age[beta.order[cutoff+1:n_ind]] <- "Low Risk - Age"

age <- factor(fac.age)
}

#cova <- data.frame(age=age, gender=clinical_info.expr[id.sele,"gender"])
cova <- data.frame(age=age)
#x11(width=10);
survit(stime=getstime(na.omit(clinical_info.expr[id.sele,]))[,1], 
event=getstime(na.omit(clinical_info.expr[id.sele,]))[,2], 
expr=na.omit(combat_data_unique.sel[gene, id.sele]), 
X=list(), 
gene.label=gene, 
data.label=t2,
boxplot.ylab="Gene Expression",
smin=smin,
covariates=cova,
plotKM=plotKM,
include.covar=include.covar
) -> expr
return(list(cnv=cnv, expr=expr))

}

sapply("RND3", function(x)survproves(id.sel=sel, t1="CNV - Untreated GBM corrected - age",t2="Expression - Untreated GBM corrected - age", gene=x, plotKM=TRUE, include.covar=TRUE)) -> tar.rnd3
sapply(tar.tgfb.sig, function(x)survproves(id.sel=sel, t1="CNV - Untreated GBM corrected - age factor",t2="Expression - Untreated GBM corrected - age factor", gene=x, plotKM=TRUE)) -> tar.tgfb.sig
sapply(tar.rhoe, function(x)survproves(id.sel=sel, t1="CNV - Untreated GBM corrected - age factor",t2="Expression - Untreated GBM corrected - age factor", gene=x)) -> tar.rhoe.f


t(apply(tar.gefs.res, 2, function(x){

expr.hr_simple <- summary(x[["expr"]][["C"]])[["conf.int"]][1]
expr.ciu_simple <- summary(x[["expr"]][["C"]])[["conf.int"]][4]
expr.cil_simple <- summary(x[["expr"]][["C"]])[["conf.int"]][3]
expr.p_simple <- summary(x[["expr"]][["C"]])[["sctest"]]["pvalue"]
cnv.hr_simple <- summary(x[["cnv"]][["C"]])[["conf.int"]][1]
cnv.ciu_simple <- summary(x[["cnv"]][["C"]])[["conf.int"]][4]
cnv.cil_simple <- summary(x[["cnv"]][["C"]])[["conf.int"]][3]
cnv.p_simple <- summary(x[["cnv"]][["C"]])[["sctest"]]["pvalue"]
expr.hr_age <- summary(x[["expr"]][["CA"]])[["conf.int"]][1,1]
expr.ciu_age <- summary(x[["expr"]][["CA"]])[["conf.int"]][1,4]
expr.cil_age <- summary(x[["expr"]][["CA"]])[["conf.int"]][1,3]
expr.p_gene <- summary(x[["expr"]][["CA"]])[["coefficients"]][1,"Pr(>|z|)"]
expr.p_model <- summary(x[["expr"]][["CA"]])[["sctest"]]["pvalue"]
cnv.hr_age <- summary(x[["cnv"]][["CA"]])[["conf.int"]][1,1]
cnv.ciu_age <- summary(x[["cnv"]][["CA"]])[["conf.int"]][1,4]
cnv.cil_age <- summary(x[["cnv"]][["CA"]])[["conf.int"]][1,3]
cnv.p_gene <- summary(x[["cnv"]][["CA"]])[["coefficients"]][1,"Pr(>|z|)"]
cnv.p_model <- summary(x[["cnv"]][["CA"]])[["sctest"]]["pvalue"]

return(c(
expr.hr=expr.hr_simple, 
expr.upper=expr.ciu_simple, 
expr.lower=expr.cil_simple, 
expr.p=expr.p_simple, 
cnv.hr=cnv.hr_simple, 
cnv.upper=cnv.ciu_simple, 
cnv.lower=cnv.cil_simple, 
cnv.p=cnv.p_simple, 
expr.hr.multi=expr.hr_age, 
expr.upper.multi=expr.ciu_age, 
expr.lower.multi=expr.cil_age, 
expr.p.gene=expr.p_gene, 
expr.p.multi=expr.p_model, 
cnv.hr.multi=cnv.hr_age, 
cnv.upper.multi=cnv.ciu_age, 
cnv.lower.multi=cnv.cil_age, 
cnv.p.gene=cnv.p_gene, 
cnv.p.multi=cnv.p_model))

})) -> tar.gefs.res.tab


