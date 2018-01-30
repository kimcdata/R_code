Bootstrap <- function(dataset1, ntargets = nrow(dataset1), dataset2=NULL, ncpu=6, ngenes=1000, method="spearman", coverage=10000){

library(snowfall)
sfInit(parallel=T, cpus=ncpu)
sfExport("method")
sfExport("ds1")
sfExport("ngenes")

ds1 <- dataset1
ds2 <- dataset2
ncor <- nrow(ds1)*ntargets
ncor <- ncor*coverage
nran <- ((ngenes*ngenes)-ngenes)/2

round(ncor/nran) -> nrep

resampled.correlation <- sfSapply(1:nrep, function(x){
dat <- ds1[sample(nrow(ds1), ngenes),]
dat <- apply(dat, 2, sample, replace=T)
cor.mat <- cor(t(dat), method=method)
cor.mat <- cor.mat[upper.tri(cor.mat)]
return(cor.mat)
})

return(resampled.correlation)

}

