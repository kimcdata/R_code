read.table("clipboard",F) -> pgc
pgc[,2] -> dat
pgc[,1] -> fac
sub("Older","Elderly",fac) -> fac
tapply(dat, fac, function(x)x) -> tmp
-as.numeric(as.factor(fac)) -> fac
boxplot(dat~fac, names=c("Young", "Elderly"))

t.test(tmp[[1]], tmp[[2]])
