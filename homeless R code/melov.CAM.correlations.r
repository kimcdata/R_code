getBootstrapDistBM <- function(data.n){
 tmp <- apply(data.n,2,sample,replace=T)
 require(bigmemory)
 bm <- filebacked.big.matrix(nrow(data.n), nrow(data.n), type="double",backingfile="data3.tmp.idx",backingpath="D:/tmp")
 res <- sapply(1:nrow(data.n), function(x) bm[x,] <<- cor(tmp[x,], t(tmp), method="spearman"))
 return(lower.tri.bm(bm))
}

lower.tri.bm <- function(x,diag=F) {
  require(biganalytics)
  col <- big.matrix(nrow(x),ncol(x),type="integer")
  row <- big.matrix(nrow(x),ncol(x),type="integer")
  tmp <- sapply(1:ncol(col),function(y) row[,y] <<- 1:nrow(col))
  tmp <- sapply(1:ncol(row),function(y) col[,y] <<- rep(y,nrow(row)))
  tf <- big.matrix(nrow(x),ncol(x),type="double")
  tmp <- sapply(1:ncol(row),function(y){
    if(diag){
      sel <- row[,y] >= col[,y]
    }else{
      sel <- row[,y] > col[,y]
    }
    tf[which(sel),y] <<- x[which(sel),y]
    tf[which(!sel),y] <<- NA
  })
  mea <- mean(as.vector(as.matrix(tf)),na.rm=T)
  sds <- sd(as.vector(as.matrix(tf)),na.rm=T)
  return(c(Mean=mea,SD=sds))
}

getBootstrapDist <- function(data.n){

 sel <- sample(1:dim(data.n)[2], dim(data.n)[2], replace=T)
 corr <- cor(t(apply(data.n,2,sample,replace=T)),method="spearman")
 #corr <- corr[lower.tri(corr)]
 #quant <- quantile(corr[lower.tri(corr)], probs=probs)
 return(c(Mean=mean(corr[lower.tri(corr)]),SD=sd(corr[lower.tri(corr)])))
 #return(corr)
}


getBootstrapBigDist <- function(data.n, size=1000){

 sel <- sample(1:dim(data.n)[1], size, replace=T)
 corr <- cor(t(apply(data.n[sel,],2,sample,replace=T)),method="spearman")
 #corr <- corr[lower.tri(corr)]
 #quant <- quantile(corr[lower.tri(corr)], probs=probs)
 return(c(Mean=mean(corr[lower.tri(corr)]),SD=sd(corr[lower.tri(corr)])))
 #return(corr)

}

bigdist <- list()
for(i in 1:10){
bigdist[[i]] <- getBootstrapBigDist(ctu.mean,size=10000)
}
ctu.mean.dist <- bigdist[[sample(1:10,1)]]

ctu.eif6.cor <- getBootStrapExistingDist(ctu.mean, ctu.eif6, ctu.mean.dist)


getBootStrapExistingDist <- function(data, variable,dist,method="spearman"){

corr <- t(cor(variable, t(data), method=method))[,1]

p.nonmulti <- c()
p.multi <- c()

for(i in 1:length(corr)){

if(corr[i]<0){

#sel <- which(dist<corr[i])
#if(length(sel)>0){
#p.nonmulti[i] <- as.numeric(unlist(strsplit(rev(names(sel))[1], "%")))/100
p.nonmulti[i] <- pnorm(corr[i], mean=dist[1], sd=dist[2])
#} else {
#p.nonmulti[i] <- 0
#}
}

if(corr[i]>0){

#sel <- which(dist>corr[i])
#if(length(sel)>0){
#p.nonmulti[i] <- 1-as.numeric(unlist(strsplit(names(sel)[1], "%")))/100
p.nonmulti[i] <- 1-pnorm(corr[i], mean=dist[1], sd=dist[2])
#} else {
#p.nonmulti[i] <- 0
#}
}

if(corr[i]==0){
p.nonmulti[i]=1
}
}

return(list(PV = p.nonmulti, correlation = corr))

}

melov.young.pnorm <- getBootstrapDist(melov.n[,which(cl=="Young")])
melov.old.pnorm <- getBootstrapDist(melov.n[,which(cl="Old")])

eif6.corr.melov.young <- getBootStrapExistingDist(melov.n[,which(cl=="Young")], melov.n[5142,which(cl=="Young")], melov.young.pnorm)

eif6.corr.melov.old <- getBootStrapExistingDist(melov.n[,which(cl=="Old")], melov.n[5142,which(cl=="Old")], melov.old.pnorm)








