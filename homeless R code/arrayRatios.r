arrayRatios <- function(data, groups, correctionMethod="BH", result="p", mu=0, rat=FALSE){

	m <- correctionMethod

	stopifnot(	
	any(tolower(result) == c("t","p","both")),
	any(tolower(m) == c("none","bh")),
	length(unique(groups))==2,
	is.factor(groups),
	is.numeric(data),
	is.matrix(data),
	length(groups)==ncol(data))

	g1 <- 1:table(groups)[1]
	g2 <- (table(groups)[1]+1):length(groups)
	comps <- as.matrix(expand.grid(g1,g2))

	ratios <- apply(comps, 1, function(x){
		return(data[,x[1]]-data[,x[2]])
	})
	if(rat){
		return(ratios)
		}

	res <- apply(ratios, 1, function(y){

		tt <- try(t.test(y, mu=mu))
		if(tolower(result)=="t"){
			return(tt$statistic)
		}
		if(tolower(result)=="p"){
			return(tt$p.value)
		}
		if(tolower(result)=="both"){
			return(c(t=tt$statistic, p=tt$p.value))
		}
	})
	
	means <- apply(ratios, 1, mean)
	
	if(result=="p" && m != "none"){
	print("adjusting p values...")
	res <- p.adjust(res, m)
	}
	if(result=="both"){
	print("adjusting p values...")
	res <- t(res)
	res[,2] <- p.adjust(res[,2],"BH")
	}
	
	return(list(res=res, means=means))


}

arrayRatios(dat.n, groups, result="both") -> tmp
writeClipboard(as.character(tmp$res[,1]))
writeClipboard(as.character(tmp$res[,2]))
writeClipboard(as.character(tmp$means))


