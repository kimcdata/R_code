removeDuplicates3 <- function(data, data.fac,tar){
	if(missing(data.fac)){
    if(tar=="G"||tolower(tar)=="genes"){
		  data.fac <- as.character(data[,2])
		  sw <- which(data.fac == "")
		  data.fac[sw] <- paste("unknown",1:length(sw),sep="")
	  	pbd <- data[,1:2]
  	  probes <- tapply(as.character(pbd[,1]),as.character(data.fac),function(x) x,simplify=F)
		  }
    if(tar=="P"||tolower(tar)=="probes"){
      data.fac <- as.character(data[,1])
    }
	}
	data.fac <- as.character(data.fac)
	probes <- rownames(data)
	probes <- tapply(as.character(probes),as.character(data.fac),function(x) x,simplify=F)
	#data <- data[,-c(1:2)]
	data <- as.matrix(data)
	data <- apply(data, 2, as.numeric)
	#data.fac <- as.character(data.fac)
	sel <- union(which(duplicated(data.fac)),which(duplicated(data.fac,fromLast=T)))
	data.m <- data[-sel,]
	rn <- as.character(data.fac)[-sel]
	data <- data[sel,]
	ctut <- by(data,data.fac[sel],function(x) {
		x <- as.matrix(x)
		if(nrow(x) > 1){
			o <- order(apply(x,1,mean),decreasing=T)[1]
			return(x[o,])
		}else{
			return(x)
		}
	})
	ctut2 <- t(sapply(ctut,function(x) x))
	data.f <- rbind(data.m,ctut2)
	rownames(data.f) <- c(rn,rownames(ctut2))
	if(tar=="G"||tolower(tar)=="genes"){
  	probes <- sapply(probes,paste,collapse=";")
  	probes <- probes[rownames(data.f)]
  	data.f <- data.frame(probes=probes,data.f)
	}
	if(tar=="P"||tolower(tar)=="probes"){
	data.f <- data.frame(data.f)
	}
	return(data.f)
}


