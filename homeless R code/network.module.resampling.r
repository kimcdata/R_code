
lapply(modules, function(x){
	mi <- match(x, modules)
	size <- tab[mi,3]
	size_y <- tab[mi, 1]
	size_o <- tab[mi, 2]
	
	ratio <- sapply(1:10000, function(y){
			res <- table(sample(edge_pool, size, replace=FALSE))
			ratio <- log2(res[1])-log2(res[2])
			return(ratio)
	})
	return(ratio)
}) -> resampling

sapply(resampling, mean) -> r.means
sapply(resampling, sd) -> r.sd

# tab = module edge counts for young and elderly

log2(tab) -> tab.log
apply(tab.log, 1, function(x)x[1]-x[2]) -> tmp # log2 ratio of young vs elderly

width<-1
barplot(tmp, width=width, col="white") -> barplot.pos

offset <- width/2

sapply(1:7, function(x){

segments(x0 = barplot.pos[x]-offset, y0 = r.means[x], x1 = barplot.pos[x]+offset, y1 = r.means[x])
segments(x0 = barplot.pos[x]-offset, y0 = r.means[x]+(r.sd[x]*2), x1 = barplot.pos[x]+offset, y1 = r.means[x]+(r.sd[x]*2), lty="dashed")
segments(x0 = barplot.pos[x]-offset, y0 = r.means[x]-(r.sd[x]*2), x1 = barplot.pos[x]+offset, y1 = r.means[x]-(r.sd[x]*2), lty="dashed")

})

plot(1:40,1:40, pch=1:40, cex=2)

boxplot(resamp.tab, ylim=c(-3.5,3.5), xaxt="none", border="grey30")
points(1:7,tmp, col="black", cex=2, pch=23, bg="white")
points(1:7,tmp, col="black", cex=2, pch=3, bg="white")
sapply(1:6, function(x)segments(x0 = x+0.5, x1=x+0.5, y0 = -5, y1=5, lty="dashed", col="grey"))
axis(1, at=1:7, labels=modules)


