#load data
read.delim("JetSet.filtered.data.Zscored.with.phys.var.symbols.txt",stringsAsFactors=F) -> jetset
jetset[,1] -> jetset.rowIds
jetset[,-c(1)] -> jetset 

by(jetset, jetset.rowIds, function(x)x[which.max(rowMeans(x)),]) -> jetset.by
t(sapply(jetset.by, function(x)x)) -> jetset.u
apply(jetset.u, 2, unlist) -> jetset.u

#fetch sample IDs for pre and post samples
colnames(jetset.u) -> sampleIds
grep("Pre", sampleIds) -> pre.i
grep("Post", sampleIds) -> post.i

#subset data to pre and post samples
jetset.u[,pre.i] -> jetset.pre
jetset.u[,post.i] -> jetset.post


#fetch participant numbers for all samples and pre and post samples
gsub("P.*","",gsub(".*_","",sampleIds)) -> sampleNos
sampleNos[pre.i] -> sampleNos.pre
sampleNos[post.i] -> sampleNos.post

#get matching participant numbers to find participants with both pre and post samples available
intersect(sampleNos.pre, sampleNos.post) -> finalSamples


#subset data for participants with both pre and post samples available
jetset.pre.f <- jetset.pre[,match(finalSamples, sampleNos.pre)]
jetset.post.f <- jetset.post[,match(finalSamples, sampleNos.post)]
gsub("P.*","",gsub(".*_","",colnames(jetset.pre.f))) -> sampleNos.pre.f
gsub("P.*","",gsub(".*_","",colnames(jetset.post.f))) -> sampleNos.post.f


#resample data i.e. generate random data
apply(jetset.post.f, 2, sample) -> jetset.post.resample
apply(jetset.pre.f, 2, sample) -> jetset.pre.resample


#create correlation matrices for data and resampled data
cor(t(jetset.pre.f), method="spearman") -> jetset.pre.spearman
cor(t(jetset.post.f), method="spearman") -> jetset.post.spearman
cor(t(jetset.post.resample),method="spearman") -> jetset.post.resample.cormat
cor(t(jetset.pre.resample),method="spearman") -> jetset.pre.resample.cormat

#create matrix of correlation differences for resampled data
jetset.pre.resample.cormat - jetset.post.resample.cormat -> jetset.resample.corrdiff

#calculate mean and sd of correlation differences
mean(jetset.resample.corrdiff[upper.tri(jetset.resample.corrdiff)]) -> corrdiff.resample.mean
sd(jetset.resample.corrdiff[upper.tri(jetset.resample.corrdiff)]) -> corrdiff.resample.sd

#create matrix of correlation differences for real data
jetset.pre.spearman - jetset.post.spearman -> corrdiff

#using pnorm, calculate the significance of each correlation difference using the distribution of correlation differences from the resampled data
apply(corrdiff, 1, function(x){
	sapply(x, function(y){
		if(y<=0){
			pnorm(y, mean=corrdiff.resample.mean, sd=corrdiff.resample.sd,lower.tail=T)
		}else{
			pnorm(y, mean=corrdiff.resample.mean, sd=corrdiff.resample.sd, lower.tail=F)
		}
	}) 
})-> corrdiff.pnorm

#load gene sets

read.delim("~/chronos/roger.kissane.angiogenesis.related.genes.txt",header=F) -> target.angio.genes
as.character(target.angio.genes[,1]) -> target.angio.genes
target.angio.genes[37] <- "THBS1"
target.angio.genes[60] <- "IGF2"
intersect(rownames(corrdiff.pnorm), target.angio.genes) -> target.angio.genes.matched
readLines("~/chronos/heritage.phys.var.ids.txt") -> phys.var.ids

#create gene set files for GSEA

apply(corrdiff.pnorm[c(phys.var.ids,target.angio.genes.matched),], 1, function(x)names(x[which(x<0.01)])) -> diffcorr.gene.sets
sapply(diffcorr.gene.sets, function(x)c(x,rep("",943-length(x)))) -> gmx
write.table(gmx, "Angio.genes.diffcorr.genesets.pv0.01.gmt.txt", sep="\t", quote=F,row.names=F)








