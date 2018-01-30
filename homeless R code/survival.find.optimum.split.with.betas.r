# age.beta <- rnd3.tab$age.at.diagnosis * tmp$coefficients[2]
# rnd3.beta <- rnd3.tab$RND3.expression * tmp$coefficients[1]

# order(rnd3.beta) -> rnd3.order
# rnd3.fac <- rep(0, nrow(rnd3.tab))
# rnd3.fac[rnd3.order[1:262]] <- 0
# rnd3.fac[rnd3.order[263:525]] <- 1

# smin=10

# sapply(1:length(rnd3.order), function(x){
	# if(x < smin){
	# return(1)
	# }
	# else if(x > length(rnd3.order)-smin){
		# return(1)
	# }
	# else {
		# rnd3.fac <- rep(0, length(rnd3.beta))
		# rnd3.fac[rnd3.order[1:x]] <- 0
		# rnd3.fac[rnd3.order[x+1:length(rnd3.order)]] <- 1
		# #summary(coxph(rnd3.surv~factor(rnd3.fac)+as.numeric(rnd3.tab$age.at.diagnosis)))$logtest[3] -> p
		# summary(coxph(rnd3.surv~factor(rnd3.fac)+age.fac))$coefficients[,"Pr(>|z|)"][1] -> p
		# return(p)
	# }
# }) -> rnd3.opti.2







################################################################################

survit <- function(stime, event, expr, X=list(), split, smin=(length(stime)/10), gene.label="", data.label="Gene Expression", boxplot.ylab="Gene Expression", covariates=c()){

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
fac[beta.order[1:cutoff]] <- "High Risk"
fac[beta.order[cutoff+1:n_ind]] <- "Low Risk"

covar <- data.frame(factor(fac), covariates)
colnames(covar)[1] <- "x"

coxph(surv.obj ~ factor(fac), method="breslow") -> cox.model
#coxph(surv.obj ~ factor(fac)*age*gender, method="breslow") -> cox.model.int
eval(parse(text=paste("cox.model.int <- coxph(surv.obj~",paste(colnames(tmp2),collapse="*"),",method=\"breslow\",data=covar)")))
survfit(surv.obj ~ factor(fac)) -> sf.model
survdiff(surv.obj ~ factor(fac)) -> sd.model

summary(cox.model.int) -> cox.summ
round(cox.summ$coefficients["xLow Risk","Pr(>|z|)"],6) -> pv
round(cox.summ$conf.int["xLow Risk","exp(coef)"],2) -> hr
round(cox.summ$conf.int["xLow Risk","lower .95"],2) -> ci.l
round(cox.summ$conf.int["xLow Risk","upper .95"],2) -> ci.u

#x11();
#par(mfrow=c(1,2))
pdf(filename=paste(gene.label,".",data.label,".pdf",sep=""))
nf <- layout(matrix(c(1,2), ncol=2), widths=c(3,2))
layout.show(nf)
plot(sf.model, col=c("red", "blue"), xlab="Days to death", ylab="Probability of survival", main=paste("Survival analysis - ",data.label,": ",gene.label,". \n HR = ",hr,". CI(L) = ",ci.l,". CI(U) = ",ci.u,". P-value = ",pv,sep=""), las=1)
legend("topright", col=c("red", "blue"), legend=c(paste("High Risk: ",cutoff,sep=""), paste("Low Risk: ", n_ind-cutoff,sep="")), pch="---")
print("High risk group")
print(cutoff)
print("Low risk group")
print(n_ind-cutoff)

expr.p <- t.test(expr[fac=="High Risk"], expr[fac=="Low Risk"])$p.value
#x11(); 
boxplot(expr~fac, main=paste("Expression: ", gene.label, "\n p-value = ",round(expr.p, 6), sep=""), las=1, ylab=boxplot.ylab)
dev.off()
return(list(C=cox.model, CA = cox.model.int, F=sf.model, D=sd.model))

}


id.selc <- intersect(id.sel, rownames(clinical_info.cnv))
id.sele <- intersect(id.sel, rownames(clinical_info.expr))

survit(stime=getstime(na.omit(clinical_info.cnv[id.selc,]))[,1], 
event=getstime(na.omit(clinical_info.cnv[id.selc,]))[,2], 
expr=na.omit(cnv_data_combat_unique.sel["RND3", id.selc]), 
X=list(), 
gene.label="RND3", 
data.label="CNV (all patients)",
boxplot.ylab="Copy Number Ratio",
) -> tmp

x11();
survit(stime=getstime(na.omit(clinical_info.expr[id.sele,]))[,1], 
event=getstime(na.omit(clinical_info.expr[id.sele,]))[,2], 
expr=na.omit(combat_data_unique.sel["RND3", id.sele]), 
X=list(), 
gene.label="RND3", 
data.label="Expression (Untreated primary GBM)",
boxplot.ylab="Gene Expression",
) -> tmp

survproves <- function(id.sel, t1, t2, smin=(length(id.sel)/10), gene){

id.selc <- intersect(id.sel, rownames(clinical_info.cnv))
id.sele <- intersect(id.sel, rownames(clinical_info.expr))
cova <- data.frame(age=clinical_info.cnv[id.selc,"age_at_initial_pathologic_diagnosis"], gender=clinical_info.cnv[id.selc,"gender"])
x11(width=10);
survit(stime=getstime(na.omit(clinical_info.cnv[id.selc,]))[,1], 
event=getstime(na.omit(clinical_info.cnv[id.selc,]))[,2], 
expr=na.omit(cnv_data_combat_unique.sel[gene, id.selc]), 
X=list(), 
gene.label=gene, 
data.label=t1,
boxplot.ylab="Copy Number Ratio",
smin=smin,
covariates=cova
) -> tmp

cova <- data.frame(age=clinical_info.expr[id.sele,"age_at_initial_pathologic_diagnosis"], gender=clinical_info.expr[id.sele,"gender"])
x11(width=10);
survit(stime=getstime(na.omit(clinical_info.expr[id.sele,]))[,1], 
event=getstime(na.omit(clinical_info.expr[id.sele,]))[,2], 
expr=na.omit(combat_data_unique.sel[gene, id.sele]), 
X=list(), 
gene.label=gene, 
data.label=t2,
boxplot.ylab="Gene Expression",
smin=smin,
covariates=cova
) -> tmp


}

sapply(tar, function(x)survproves(id.sel=sel, t1="CNV - Untreated GBM corrected",t2="Expression - Untreated GBM corrected", gene=x))
#survit(stime, status, rnd3.tab$RND3.expression, X=list(age=age.fac2), gene.label="RND3") -> tmp


#### integrating additional information into the cox model

### RAC1 example gene


id.sele <- intersect(id.sel, rownames(clinical_info.expr))
gene <- "RND3"
surv.obj <- Surv(getstime(na.omit(clinical_info.expr[id.sele,]))[,1], getstime(na.omit(clinical_info.expr[id.sele,]))[,2])
rac1.beta <- coxph(surv.obj~na.omit(combat_data_unique.sel[gene, id.sele]))$coefficients[1]
rac1.beta.order <- order(rac1.beta*combat_data_unique.sel[gene, id.sele], decreasing=T)
n_ind <- length(rac1.beta.order)

sapply(1:n_ind, function(x){
	if(x < smin){
	return(1)
	}
	else if(x > n_ind-smin){
		return(1)
	}
	else {
		fac <- rep(0, n_ind)
		fac[rac1.beta.order[1:x]] <- 0
		fac[rac1.beta.order[x+1:n_ind]] <- 1
		#summary(coxph(rnd3.surv~factor(rnd3.fac)+as.numeric(rnd3.tab$age.at.diagnosis)))$logtest[3] -> p
		summary(coxph(surv.obj~factor(fac), method="breslow"))$coefficients[,"Pr(>|z|)"][1] -> p
		return(p)
	}
}) -> p.opti

which.min(p.opti) -> cutoff

fac <- rep(0, n_ind)
fac[rac1.beta.order[1:cutoff]] <- "High Risk"
fac[rac1.beta.order[cutoff+1:n_ind]] <- "Low Risk"

age <- clinical_info.expr[id.sele,"age_at_initial_pathologic_diagnosis"]
gender <- factor(clinical_info.expr[id.sele,"gender"])
kps <- as.numeric(clinical_info.expr[id.sele,"karnofsky_performance_score"])

coxph(surv.obj ~ factor(fac)*age*gender, method="breslow") -> cox.model.int



















