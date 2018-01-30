############# Affymetrix microarray preprocessing ###################

library(affy)
library(simpleaffy)
library(gcrma)
library(preprocessCore)
library(gplots)
library(biomaRt)

files <- dir("cel/",full.names = T)

raw <- read.affybatch(files)
detection.p.val(raw) -> raw.calls
rowSums(raw.calls$call == "P") -> raw.pcount

names(raw.pcount[raw.pcount > round(length(colnames(dat.rma))*0.5,0)]) -> probes.sel

exprs(rma(raw)) -> dat.rma
gsub("_.*","",colnames(dat.rma)) -> colnames(dat.rma)

dat.rma[probes.sel,] -> dat.filt

#affy2symbol(rownames(dat.filt)) -> probe.symbols
#probe.symbols$hgnc_symbol[match(rownames(dat.filt), probe.symbols$affy_hg_u133_plus_2)] -> dat.hgnc

read.delim("GPL570-55999.txt",header = T,sep = "\t",stringsAsFactors = FALSE,skip=16) -> array.annot

array.annot$Gene.Symbol[match(rownames(dat.filt), array.annot$ID)] -> dat.hgnc

by(dat.filt, dat.hgnc, function(x)x[which.max(rowSums(x)),]) -> dat.by
t(sapply(dat.by, function(x)x)) -> dat.u
apply(dat.u, 2, unlist) -> dat.u
dat.u <- dat.u[-1,]

read.delim("GSE25507_sample_info.txt",stringsAsFactors = FALSE) -> sample.info


