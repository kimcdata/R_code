overlapper <- function(exp1, exp1.name="A", exp2, exp2.name="B", exp1.full, exp2.full, filename=paste("venn.diagram.",unlist(strsplit(as.character(Sys.time()), " "))[1],".tif",sep="")){

library(VennDiagram)
exp1.full.unique <- toupper(unique(exp1.full))
exp2.full.unique <- toupper(unique(exp2.full))
full.inter <- intersect(exp1.full.unique, exp2.full.unique)
exp1.m <- toupper(unique(exp1))
exp2.m <- toupper(unique(exp2))
exp1.m <- full.inter[na.omit(match(exp1.m, full.inter))]
exp2.m <- full.inter[na.omit(match(exp2.m, full.inter))]
venn.diagram(list(exp1.name=exp1.m, exp2.name=exp2.m), sub=paste("p-value: ",getP(res=length(intersect(exp1.m, exp2.m)), length(exp1.m), length(exp2.m), length(full.inter))$PN,sep=""), filename=filename)
return(list(exp1.name= exp1.m[!exp1.m %in% exp2.m], shared=exp1.m[exp1.m %in% exp2.m], exp2.name=exp2.m[!exp2.m %in% exp1.m]))
}


overlapper2 <- function(exp1, exp2,exp1.full, exp2.full, lower.tail=TRUE){

exp1.full.unique <- toupper(unique(exp1.full))
exp2.full.unique <- toupper(unique(exp2.full))
full.inter <- intersect(exp1.full.unique, exp2.full.unique)
exp1.m <- toupper(unique(exp1))
exp2.m <- toupper(unique(exp2))
exp1.m <- full.inter[na.omit(match(exp1.m, full.inter))]
exp2.m <- full.inter[na.omit(match(exp2.m, full.inter))]

x <- length(intersect(exp1.m, exp2.m))
m <- length(exp1.m)
n <- length(full.inter)-length(exp1.m)
k <- length(exp2.m)
return(list(int = length(intersect(exp1.m, exp2.m)),stat = phyper(q=x,m=m,n=n,k=k, lower.tail=lower.tail), genes = intersect(exp1.m, exp2.m)))
}

eif6.up.hypox.down.overlapper <- overlapper(exp1=eif6.up, exp1.name="EIF6 +/- Up-regulated", exp2=hypox.down, exp2.name="Hypoxia Down-regulated", filename="EIF6.up.hypox.down.venn.tiff", exp1.full = eif6.full, exp2.full=hypox.full)

eif6.up.hypox.up.overlapper <- overlapper(exp1=eif6.up, exp1.name="EIF6 +/- Up-regulated", exp2=hypox.up, exp2.name="Hypoxia Up-regulated", filename="EIF6.up.hypox.up.venn.tiff", exp1.full = eif6.full, exp2.full=hypox.full)

eif6.down.hypox.down.overlapper <- overlapper(exp1=eif6.down, exp1.name="EIF6 +/- Down-regulated", exp2=hypox.down, exp2.name="Hypoxia Down-regulated", filename="EIF6.down.hypox.down.venn.tiff", exp1.full = eif6.full, exp2.full=hypox.full)

eif6.down.hypox.up.overlapper <-overlapper(exp1=eif6.down, exp1.name="EIF6 +/- Down-regulated", exp2=hypox.up, exp2.name="Hypoxia Up-regulated", filename="EIF6.down.hypox.up.venn.tiff", exp1.full = eif6.full, exp2.full=hypox.full)


getP <- function(res,x,y, max,rep=100000){

	t1 <- c()

	for(i in 1:rep){
		t1[i] <- length(intersect(sample(max,x, replace=F), sample(max,y, replace=F)))
	}

	f <- which(quantile(t1, seq(0,1,0.0001))>res)
	if(length(f)==0){
		p.quant <- 0
	} else {
		p.quant <- 1-(as.numeric(sub("%","",names(f)[1]))/100)
	}

	#p.dist <- shapiro.test(t1)
	pval <- pnorm(res, mean=mean(t1), sd=sd(t1),lower.tail=FALSE)
	print(paste("mean=",mean(t1)))
	print(paste("pval=",p.quant))

	return(list(PN = pval, PQ = p.quant))
}


overlapper.test <- lapply(1:1000, function(x){

sample(10000, x) -> t1
sample(10000, x) -> t2
return(overlapper2(t1,t2,1:3000, 1:3000, lower.tail=FALSE)$stat)

})