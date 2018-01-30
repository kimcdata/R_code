
# resampling using snowfall
# creates correlation matrices from random subsets of 1000 genes at a time

library(snowfall)
ng <- 1000 #number of genes
ncpu <- 6 # number of cpus
nrep <- 1000 # number of random subsets
sfInit(parallel=T, cpus=ncpu)
sfExport("copd.data") #export your data to the snowfall processes

random.correlation.values <- sfLapply(1:nrep, function(x){

copd.rand <- apply(copd.data, 2, sample, replace=T) # randomise each column in the data

row.sel <- sample(nrow(copd.rand), 1000, replace=F) # select 1000 random rows

cor.mat <- cor(t(copd.rand[row.sel,]), method="spearman") # calculate correlation matrix using 1000x1000 subset

cor.mat <- cor.mat[upper.tri(cor.mat)] # get upper tri

return(as.vector(cor.mat))

})


quantile.distribution <- quantile(random.correlation.values, probs=seq(0,1,0.0001))

# this code matches the real correlation values to the quantiles of the random values, finds it's position and uses the numeric "name" of that position i.e. "99.7737%" = 0.997737 to calculate the p value. e.g 1-0.997737 = 0.002.... 
# run the command "names" on the results of the quantile function and you will see what i mean.
# the sub is to remove the % sign, the /100 is to convert the % into a value between 0 and 1
# the code for correlation values > 0 is different to the code for < 0. For > 1 you need to match the first quantile that is above your correlation value, hence names(f)[1]. For < 0 you need to match the last quantile value that is less than your correlation value, hence names(f)[length(whcih(qu<x))]

quantile.p.value <- sapply(real.cor.values, function(x){
	if(x>0){
		f <- which(quantile.distribution>x)
		if(length(f)==0){
			p.quant <- 0
		} else {
			p.quant <- 1-(as.numeric(sub("%","",names(f)[1]))/100)
		}
	} else {
		f <- which(quantile.distribution<x)
		if(length(f)==0){
			p.quant <- 0
		} else {
			p.quant <- as.numeric(sub("%","",names(f)[length(which(qu<x))]))/100
		}
	}
	return(p.quant)
})

