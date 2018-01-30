HOPACH2LIST <- function(data, dmat=NULL, d="cor", clusters = "best", K=15, kmax=9, khigh=9, coll="seq", newmed="medsill",mss = "med",impr = 0, initord = "co", ord = "own", verbose=FALSE, cluster.size.cutoff = 10){

hopach(data=data, dmat=dmat, d=d, clusters = clusters, K=K, kmax=kmax, khigh=khigh, coll=coll, newmed=newmed,mss = mss,impr = impr, initord = initord, ord = ord, verbose=verbose) -> data.hopach

clusters <- tapply(1:nrow(data), data.hopach$clustering$labels, function(x)x)
sapply(clusters, length) -> clust.sizes
clusters.cutoff <- names(which(clust.sizes>=cluster.size.cutoff))
clusters.cutoff.data <- clusters[clusters.cutoff]

res <- lapply(names(clusters.cutoff.data), function(x){
expr <- data[clusters.cutoff.data[[x]],]
med <- data[data.hopach$clustering$medoids[na.omit(match(clusters.cutoff.data[[x]],data.hopach$clustering$medoids))],]
return(list(expression=expr, medoid=med))
})
return(res)

}

sapply(tum.greedy, function(x){

x11()
dat <- x$expression

heatmap.2(dat, scale="row", col=greenred(50), Colv=FALSE, Rowv=TRUE, dendrogram="row", trace="none", labRow="", distfun = function(dat)as.dist(1-cor(t(dat))), hclustfun = function(dat)hclust(dat, "average"), margins=c(15,5), main=nrow(dat))


})

dat <- data.sig[match(m.down, sig.genes),]
heatmap.2(dat, scale="row", col=greenred(50), Colv=FALSE, Rowv=TRUE, dendrogram="row", trace="none", labRow="", distfun = function(dat)as.dist(1-cor(t(dat))), hclustfun = function(dat)hclust(dat, "average"), margins=c(15,5), main=paste("Expression in H2O2 treated C2C12s of\n ",nrow(dat)," EIF6 targets down-regulated\n in the muscle",sep=""))


heatmap.2(tum.4, scale="row", col=greenred(50), Colv=FALSE, Rowv=TRUE, dendrogram="row", trace="none", distfun = function(dat)as.dist(1-cor(t(dat))), hclustfun = function(dat)hclust(dat, "average"), main=nrow(tum.4))


t(Make.Z(t(tmp))) -> tmp
plot(1:109, tmp[4,], type="l", lwd=2, ylim=c(min(tmp), max(tmp)))
lines(1:109, tmp[2,], type="l", lwd=2)
lines(1:109, tmp[3,], type="l", lwd=2)
lines(1:109, tmp[5,], type="l", lwd=2)
