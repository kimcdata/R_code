HOPACH2LIST <- function(data, dmat=NULL, d="cor", clusters = "best", K=15, kmax=9, khigh=9, coll="seq", newmed="medsill",mss = "med",impr = 0, initord = "co", ord = "own", verbose=FALSE, cluster.size.cutoff = 10){

hopach(data=data, dmat=dmat, d=d, clusters = clusters, K=K, kmax=kmax, khigh=khigh, coll=coll, newmed=newmed,mss = mss,impr = impr, initord = initord, ord = ord, verbose=verbose) -> data.hopach

clusters <- tapply(1:nrow(data), data.hopach$clustering$labels, function(x)x)
sapply(clusters, length) -> clust.sizes
clusters.cutoff <- names(which(clust.sizes>cluster.size.cutoff))
clusters.cutoff.data <- clusters[clusters.cutoff]

res <- lapply(names(clusters.cutoff.data), function(x){
expr <- data[clusters.cutoff.data[[x]],]
med <- data[data.hopach$clustering$medoids[na.omit(match(clusters.cutoff.data[[x]],data.hopach$clustering$medoids))],]
return(list(expression=expr, medoid=med))
})
return(res)

}

sapply(data.h2l, function(x){

x11()
dat <- x$expression

heatmap.2(dat, scale="row", col=greenred(50), Colv=FALSE, Rowv=TRUE, dendrogram="row", trace="none", labRow=lab, distfun = function(dat)as.dist(1-cor(t(dat))), hclustfun = function(dat)hclust(dat, "average"), margins=c(15,5), main=nrow(dat), key=FALSE)


})