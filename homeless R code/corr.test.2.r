# int = integrated
# r.e = raue exercise
# r.a = raue ageing
# stone = stone muscle ageing
# m = melov ageing + exercise

getCorrTest <- function(ndat, gn, tar, dup=1){

sel <- which(tolower(gn)==tolower(tar))
dat <- ndat

if(length(sel)>1){
	sel <- sel[dup]
	cor.p <- apply(dat, 1, function(x){
		return(cor.test(dat[sel,], x, type="spearman")$p.value)
	})
	cor.r <- apply(dat, 1, function(x){
		return(cor.test(dat[sel,], x, type="spearman")$estimate)
	})
	} else {
		cor.p <- apply(dat, 1, function(x){
			return(cor.test(dat[sel,], x, type="spearman")$p.value)
		})
		cor.r <- apply(dat, 1, function(x){
			return(cor.test(dat[sel,], x, type="spearman")$estimate)
		})
}
return(list(P=p.adjust(cor.p, "BH"), R=cor.r))



}

getCorrTest(norm.filt[,o.c1],gn,"EIF6") -> eif6.o.c1.corr
getCorrTest(norm.filt[,o.c2],gn,"EIF6") -> eif6.o.c2.corr
getCorrTest(norm.filt[,y.c1],gn,"EIF6") -> eif6.y.c1.corr
getCorrTest(norm.filt[,y.c2],gn,"EIF6") -> eif6.y.c2.corr

getCorrTest(norm.filt[,o.c1],gn,"EIF3K",1) -> eif3k.o.c1.corr
getCorrTest(norm.filt[,o.c2],gn,"EIF3K",1) -> eif3k.o.c2.corr
getCorrTest(norm.filt[,y.c1],gn,"EIF3K",1) -> eif3k.y.c1.corr
getCorrTest(norm.filt[,y.c2],gn,"EIF3K",1) -> eif3k.y.c2.corr
 
