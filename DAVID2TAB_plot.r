# examples 

DAVID2TABLE("sol.up.human.david.txt", len=sol.up.length, tar=c("GOTERM_BP_FAT","KEGG_PATHWAY")) -> tab1
DAVID2TABLE("human.corr.neg.david.txt", len=sol.up.length, tar=c("GOTERM_BP_FAT","KEGG_PATHWAY")) -> tab2
david.comp <- COMBINETABLE(tab1, tab2, sig.colour="red", nonsig.colour="grey", force.scale=FALSE, thresh=0.1, tabx.pheno="sol.up", taby.pheno="human.corr.neg", xlab="Soleus - up regulated", ylab="Negative correlation in human model", write.tab=TRUE,write.filename="sol.up.human.corr.neg.2wayPathway.tab.txt")


DAVID2TABLE("sol.down.human.david.txt", len=sol.up.length, tar=c("GOTERM_BP_FAT","KEGG_PATHWAY")) -> tab1
DAVID2TABLE("human.corr.pos.david.txt", len=sol.up.length, tar=c("GOTERM_BP_FAT","KEGG_PATHWAY")) -> tab2
david.comp <- COMBINETABLE(tab1, tab2, sig.colour="red", nonsig.colour="grey", force.scale=TRUE, thresh=0.1,tabx.pheno="sol.down", taby.pheno="human.corr.pos", xlab="Soleus - down regulated", ylab="Positive correlation in human model", write.filename="sol.down.human.corr.pos.2wayPathway.tab.txt", write.tab=FALSE)

# functions

COMBINETABLE <- function(tabx, taby,thresh=1, sig.colour = "red", nonsig.colour = "black", filename="two.way.david.image.tif", height=480, width=480, compression="none",main="2-way Pathway Enrichment", xlab="", ylab="", legpos="topright", force.scale=TRUE, tabx.pheno = "DAVID file 1", taby.pheno= "DAVID file 2", write.tab=FALSE, write.filename="TwoWayDavidTable.txt", write.tif = FALSE, sig.choice="any"){

	repl <- c(0,1)
	names(repl) <- c("count.norm","benj")
	terms <- unique(c(tabx[,1], taby[,1]))
	
	tabx[match(terms, tabx[,1]),][,4] -> tabx.c
	taby[match(terms, taby[,1]),][,4] -> taby.c
	
	tabx[match(terms, tabx[,1]),][,2:3] -> tabx
	taby[match(terms, taby[,1]),][,2:3] -> taby
	
	rownames(tabx) <- rownames(taby) <- terms
	cbind(tabx, taby) -> tabf
	tabf.na <- which(is.na(tabf),T)
	tabf.na <- tapply(tabf.na[,1],tabf.na[,2],function(x) x)
	for(i in 1:length(tabf.na)) tabf[tabf.na[[i]],as.numeric(names(tabf.na)[i])] <- repl[colnames(tabf)[as.numeric(names(tabf.na)[i])]]
	tabf.n <- apply(tabf, 2, as.numeric)
	rownames(tabf.n) <- terms
	if(sig.choice=="any"){
		apply(tabf.n, 1, function(x)min(c(x[2], x[4]))) -> sig
	} else if(sig.choice=="both"){
		apply(tabf.n, 1, function(x){
			return(ifelse(length(which(c(x[2], x[4])<thresh))>1, sig.colour, nonsig.colour))
		})
	}
	sapply(sig, function(x)if(x<thresh){return(sig.colour)}else{ return(nonsig.colour)}) -> point.colours
	which(table(which(tabf.n[,c(1,3)]>0,T)[,1])>1) -> mat
	if(write.tif){
		tiff(filename, height=height, width=width, compression=compression)
	}
	if(force.scale){
		xlim <- ylim <- max(c(tabf.n[tabf.n[,1]<1,1], tabf.n[tabf.n[,3]<1,3]))
	} else {
		xlim <- max(tabf.n[tabf.n[,1]<1,1])
		ylim <- max(tabf.n[tabf.n[,3]<1,3])
	}
	
	#x11();plot(tabf.n[,1], tabf.n[,3], col=point.colours, pch=19, xlim=c(0,xlim), ylim=c(0,ylim), main=main, xlab=xlab, ylab=ylab)
	plot(tabf.n[,1], tabf.n[,3], col=point.colours, pch=19, xlim=c(0,xlim), ylim=c(0,ylim), main=main, xlab=xlab, ylab=ylab)
	tmp.lm <- lm(tabf.n[mat,3]~tabf.n[mat,1])
	pv <- summary.lm(tmp.lm)$coefficients["tabf.n[mat, 1]","Pr(>|t|)"]
	abline(tmp.lm, lty="dashed", col="darkgrey")
	legend(legpos, legend=c(paste("p-value: ",signif(summary.lm(tmp.lm)$coefficients["tabf.n[mat, 1]","Pr(>|t|)"],4),sep=""), paste("r squared: ",round(summary.lm(tmp.lm)$r.squared,3),sep="")),y.intersp=1.5)
	if(write.tif){
		dev.off()
	}
	tabx.c[is.na(tabx.c)] <- 0
	taby.c[is.na(taby.c)] <- 0
	tabf.n <- cbind(tabf.n[,1],tabx.c, tabf.n[,2:3], taby.c, tabf.n[,4])
	apply(tabf.n, 2 ,as.numeric) -> tabf.n
	colnames(tabf.n) <- c(paste(tabx.pheno, "normalised count"),paste(tabx.pheno, "count"),paste(tabx.pheno, "benjamini"),paste(taby.pheno, "normalised count"),paste(taby.pheno, "count"),paste(taby.pheno, "benjamini"))
	rownames(tabf.n) <- terms
	intersect(which(tabf.n[,3]<1), which(tabf.n[,6]<1)) -> s1
	setdiff(which(tabf.n[,3]<1), s1) -> s2
	setdiff(which(tabf.n[,6]<1), s1) -> s3
	tabf.n[c(s1,s2,s3),] -> tabf.n
	
	if(write.tab){
		write.table(tabf.n, write.filename, sep="\t", quote=FALSE)
	}
	return(tabf.n)
}

DAVID2TABLE <- function(file, tar=c("GOTERM_BP_FAT", "KEGG_PATHWAY"), len=1, statistic="benj"){

tab <- readLines(file)
category <- sapply(tar, function(x){
return(tab[grep(x, tab)])
})

res <- lapply(category, function(x){
term <- sapply(strsplit(x, "\t"), function(x)x[2])
count <- as.numeric(sapply(strsplit(x, "\t"), function(x)x[3]))
count.norm <- as.numeric(sapply(strsplit(x, "\t"), function(x)x[3]))/len
benj <- sapply(strsplit(x, "\t"), function(x)x[switch(statistic, "benj"=12, "fdr"=13, "p"=11)])
return(list(term=term, count.norm = as.numeric(count.norm), benj=as.numeric(benj),count=as.numeric(count)))
})

if(length(tar)>1){
do.call(rbind,sapply(res, function(x)do.call(cbind,x))) -> res.tab
} else {
do.call(rbind, res) -> res.tab
t(apply(res.tab, 1, unlist)) -> res.tab
}
return(res.tab)
}