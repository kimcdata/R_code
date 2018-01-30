read.adj <- function(file,skip=T){
	ch <- readLines(file)
	if(skip){
		sel <- grep(">",ch)
		ch <- ch[-sel]
	}
	d1 <- sapply(ch,strsplit,"\t")
	#d2 <- sapply(d1,function(x) x)
	return(d1)
}

write.adj <- function(fh1, fh2){

  sif <- adj.to.sif(1, fh1)
  sx <-  correctsif(sif, fh2)
  return(sx)
}

adj.to.sif <- function(x, fh=""){
  if(grep("adj",fh)==1){
    x <- read.adj(fh, skip=T)
  }
	mm <- adj.to.mat(x)
	sif <- mat.to.sif(mm)
	return(sif)
}

correctsif <- function(x, fh){

x[-which(x[,3]==0),] -> sx

sx[,1] -> from
sx[,2] -> to

genes[match(from, probes)] -> from.genes
which(from.genes=="---") -> sel
from.genes[sel] <- from[sel]

genes[match(to, probes)] -> to.genes
which(to.genes=="---") -> sel
to.genes[sel] <- to[sel]

sx[,1] <- from.genes
sx[,2] <- to.genes

write.sif(sx, fh)

return(sx)
}

fullmat <- function(fh, i, probes, s, res){

fh <- "old_output"
i <- 12

for(i in 1:12) {

mat <<- matrix(0, nrow=length(s[,i]), ncol=length(probes))
colnames(mat) <<- probes
rownames(mat) <<- probes[s[,i]]

x <- read.adj(paste(fh,i,".adj",sep=""))

sapply(x, rod.p)

write.table(mat, "fullmat.output.test.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

}  
#return(mat)
}

rod.p <- function(y){
		y1 <- y[1]
		y <- y[-1]
		y2 <- na.tt(y[seq(1,length(y),2)])
		y3 <- na.rm(as.numeric(y[seq(2,length(y),2)]))
		mat[y1,y2] <<- y3   
		#mat[y2,y1] <<- y3
}

adj.to.mat <- function(x){
	na.rm <- function(x) x[!is.na(x)]
	na.tt <- function(x) x[which(x != "")]
	if(is.list(x)){
#	  gnames <- as.character(sapply(x,function(x) x[1]))
	  gnames <- unique(unlist(sapply(x,function(x) x[c(1,seq(2,length(x)-1,2))])))
  }else{
    gnames <- x[,1]
  }
	mm <- matrix(0,nrow=length(gnames),ncol=length(gnames))
	rownames(mm) <- colnames(mm) <- gnames
	rod <- function(y){
		y1 <- y[1]
		y <- y[-1]
		y2 <- na.tt(y[seq(1,length(y),2)])
		y3 <- na.rm(as.numeric(y[seq(2,length(y),2)]))
		mat[y1,y2] <<- y3   
		mat[y2,y1] <<- y3
	}
	if(is.list(x)){
    sapply(x,rod)
	}else{
    apply(x,1,rod)
  }
	return(mm)
}

mat.to.sif <- function(x,diag=T){
  if(diag){
	 tms <- (length(rownames(x))):1
	 a <- rep(colnames(x),tms)
	 b <- unlist(sapply(1:nrow(x),function(y) rownames(x)[y:nrow(x)]))
	 d <- as.numeric(x[lower.tri(x,diag=T)])
	}else{
	tms <- (length(rownames(x))-1):1
	 a <- rep(colnames(x)[-length(colnames(x))],tms)
	 b <- unlist(sapply(2:nrow(x),function(y) rownames(x)[y:nrow(x)]))
	 d <- as.numeric(as.dist(x))
	}
	cyto <- cbind(b,a,d)
	return(cyto)
}

write.sif <- function(x,file,...){
	write.table(x,file=file,sep="\t",quote=F,row.names=F,col.names=F,...)
}
