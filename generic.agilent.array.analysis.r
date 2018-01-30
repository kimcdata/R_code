################ Agilent array analysis #################

setRepositories(ind=1:10)
library(preprocessCore)
library(affy)
library(marray)
library(RColorBrewer)
library(AnnotationDbi)
source("https://bioconductor.org/biocLite.R")
biocLite("hgug4112a.db")
library(hgug4112a.db)
library(sva)



############### Read data ##########################

# get files paths
afd <- "agilent.files" # change this #
files <- file.path(afd,dir(afd))

# read raw data (one colour)
read.Agilent(files, name.Rf = NULL, name.Rb = NULL) -> raw

# extract foreground signal and probe labels
raw@maGf -> raw.expr
raw@maGnames@maLabels -> probes
rownames(raw.expr) <- probes

# read array annotation file
aaf <- "GPL13912-20417.txt" # change this #
read.delim(aaf,stringsAsFactors = FALSE,comment.char = "#") -> annot
# read sample annotation file
saf <- "sample.info.txt" # change this #
read.delim(saf, sep = "\t", stringsAsFactors = F) -> sample.info

# check the strucuture of each dataset

str(raw.expr)
str(annot)
str(sample.info)


################ Set up some sample order and colour palettes ###############

order(sample.info$Day,sample.info$SampleName) -> ord


############### Raw data - boxplots ###################

x11()
boxplot(log2(raw.expr[,]), las=2, names=sample.info$SampleName)


################# Splitting data and sample info dfs into Tumour and Cell Culture Samples ###########################

raw.expr[,ord] -> raw.expr
sample.info[ord,] -> sample.info

################### Setting up some expression and variance filters ########################

#remove control probes and log2 transform

expr.probes <- grep("^A_",probes)
log2(raw.expr[expr.probes,]) -> raw.expr.nc #[n]o [c]ontrols

#check background signal 
boxplot(log2(raw@maGb))

#select probes for which the median signal is well above background
apply(raw.expr.nc, 1, function(x)median(x)) -> exprfilt
which(exprfilt > 6.5) -> expr.sel

apply(raw.expr.nc, 1, sd) -> expr.sd
which(expr.sd > quantile(expr.sd)[2]) -> sd.sel

intersect(expr.sel, sd.sel) -> probe.selection

raw.expr.nc[probe.selection,] -> raw.expr.nc.filt

boxplot(raw.expr.nc.filt, las=2)

############### Quantile normalisation ##################

dimnames(raw.expr.nc.filt) -> m
normalize.quantiles(raw.expr.nc.filt) -> data.quant.norm
dimnames(raw.expr.nc.filt) <- m

boxplot(data.quant.norm, las=2)

##################### Principle Component Analysis ########################

x11()
prcomp(t(data.quant.norm)) -> pca
plot(pca$x, pch=19, cex=2)
text(x = pca$x[,1]+8, y=pca$x[,2]+2,sample.info$SampleName)


######################## Summarising by Gene Symbol and write data to file ############################

grep("^A_",rownames(vivo.quant2)) -> sel
vivo.quant2[sel,] -> vivo.quant2

annot$GeneSymbol[match(rownames(vivo.quant2), annot$ProbeID)] -> vivo.quant2.genes
by(vivo.quant2, vivo.quant2.genes, function(x)x[which.max(rowMeans(x)),]) -> vivo.by
t(sapply(vivo.by, function(x)x)) -> vivo.u
apply(vivo.u, 2, unlist) -> vivo.u
str(vivo.u)

str(vivo.sample.filt.r2)
match(colnames(vivo.u), vivo.sample.filt.r2$ColumnName)

write.table(vivo.u, file="Data.quantile.filtered.inVivo.txt", sep="\t", quote=FALSE, col.names=NA)
write.table(t(vivo.sample.filt.r2), file="SampleInfo.quantile.filtered.inVivo.txt", sep="\t", quote=FALSE, col.names=NA)


