
grep("^A_", raw@maGnames@maLabels) -> sel
new("RGList", list(R=raw@maRf[sel,], G=raw@maGf[sel,], Rb=raw@maRb[sel,], Gb=raw@maGb[sel,], genes=raw@maGnames@maLabels[sel])) -> torin.rg.obj


data.loess.nobc
data.loess.bc
data.loess.nobc.quant
data.loess.bc.quant


normalizeWithinArrays(prol.rg.obj, method="loess", bc.method="none") -> prol.data.loess.nobc
normalizeWithinArrays(diff.rg.obj, method="loess", bc.method="none") -> diff.data.loess.nobc
normalizeWithinArrays(torin.rg.obj, method="loess", bc.method="none") -> torin.data.loess.nobc

# normalizeWithinArrays(prol.rg.obj, method="loess", bc.method="subtract") -> prol.data.loess.nobc
# normalizeWithinArrays(diff.rg.obj, method="loess", bc.method="subtract") -> diff.data.loess.nobc
# normalizeWithinArrays(torin.rg.obj, method="loess", bc.method="subtract") -> torin.data.loess.nobc

normalizeWithinArrays(prol.rg.obj, method="none", bc.method="none") -> prol.data.loess.nobc
normalizeWithinArrays(diff.rg.obj, method="none", bc.method="none") -> diff.data.loess.nobc
normalizeWithinArrays(torin.rg.obj, method="none", bc.method="none") -> torin.data.loess.nobc

# normalizeWithinArrays(prol.rg.obj, method="none", bc.method="subtract") -> prol.data.loess.nobc
# normalizeWithinArrays(diff.rg.obj, method="none", bc.method="subtract") -> diff.data.loess.nobc
# normalizeWithinArrays(torin.rg.obj, method="none", bc.method="subtract") -> torin.data.loess.nobc

par(mfrow=c(1,3))
sapply(1:3, function(x){

plotMA(prol.data.loess.nobc, array=x, status=1)
abline(h=0, col="red", lty="dashed")

})

rownames(prol.data.loess.nobc$M) <- rownames(diff.data.loess.nobc$M) <- rownames(torin.data.loess.nobc$M) <- prol.data.loess.nobc$genes
 
apply(prol.data.loess.nobc$A, 1, mean) -> loess.nobc.means
which(loess.nobc.means>6.5) -> sel
write.table(prol.data.loess.nobc$M[sel,], file="prol.data.loess.nobc.6.5.txt", sep="\t", quote=FALSE)

apply(diff.data.loess.nobc$A, 1, mean) -> loess.nobc.means
which(loess.nobc.means>6.5) -> sel
write.table(diff.data.loess.nobc$M[sel,], file="diff.data.loess.nobc.6.5.txt", sep="\t", quote=FALSE)

apply(torin.data.loess.nobc$A, 1, mean) -> loess.nobc.means
which(loess.nobc.means>6.5) -> sel
write.table(torin.data.loess.nobc$M[sel,], file="torin.data.loess.nobc.6.5.txt", sep="\t", quote=FALSE)

normalizeWithinArrays(rg.obj, method="loess", bc.method="subtract") -> data.loess.bc
normalizeWithinArrays(rg.obj, method="none", bc.method="none") -> data.malist.nobc
normalizeWithinArrays(rg.obj, method="none", bc.method="subtract") -> data.malist.bc
normalizeBetweenArrays(data.loess.nobc, method="Aquantile") -> data.loess.nobc.quant
normalizeBetweenArrays(data.loess.bc, method="Aquantile") -> data.loess.bc.quant

rownames(data.loess.nobc$M) <- rownames(data.loess.bc$M) <- rownames(data.loess.nobc.quant$M) <- rownames(data.loess.bc.quant$M) <- data.loess.nobc$genes

cols <- c("15wt vs 18het","34wt vs 19het","12wt vs 19het","13wt vs 16het","4wt vs 3het","5wt vs 4het","6wt vs 5het","3wt vs 4het")

colnames(data.loess.nobc$M) <- colnames(data.loess.bc$M) <- colnames(data.loess.nobc.quant$M) <- colnames(data.loess.bc.quant$M) <- colnames(data.malist.bc$M) <- colnames(data.malist.nobc$M) <- cols

apply(data.loess.nobc$A, 1, mean) -> loess.nobc.means
which(loess.nobc.means>7) -> sel
write.table(data.loess.nobc$M[sel,], file="data.loess.nobc.7.txt", sep="\t", quote=FALSE)

write.table(data.loess.nobc$M, file="data.loess.nobc.txt", sep="\t", quote=FALSE)
write.table(data.loess.bc$M, file="data.loess.bc.txt", sep="\t", quote=FALSE)
write.table(data.loess.nobc.quant$M, file="data.loess.nobc.quant.txt", sep="\t", quote=FALSE)
write.table(data.loess.bc.quant$M, file="data.loess.bc.quant.txt", sep="\t", quote=FALSE)



