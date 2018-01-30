source("Y:\\Normal Muscle Model\\getCels.r")

test <- david_assign2(conversion, big_data.combat)
young.old.unique2 <- cbind(Genes = rownames(young.old.unique), Probes=as.vector(unlist(attributes(young.old.unique)[3])), young.old.unique)
cbind(genes[match(old.unique.writeable[,2], probes)],old.unique.writeable) -> old.unique.writeable
write.table(old.unique.writeable, file="U133A.Plus2.mas5.OLD.unique.aracne.txt", sep="\t", quote=FALSE, row.names=FALSE)


test2 <- test2[-grep("AFFX",test2[,2]),]
  

prcomp2 <- function(input, tar, pos="bottomleft"){
sapply(as.numeric(as.factor(tar)), function(x)rainbow(length(levels(tar)))[x]) -> cols
prcomp(t(input),scale=T) -> pca
plot(pca$x, col=cols, pch=19)
legend(pos, as.character(unique(tar)), col=unique(cols), pch=19)
return(pca)
}

getCels <- function(cel.list, dataset.list=c(), norm=TRUE, platform=FALSE){

library(affy)
load("Y:/Normal Muscle Model/Raw files/annotation/conversion.RData")
source("Y:\\Normal Muscle Model\\Raw files\\annotation\\david_assign2.r")


files <- sapply(dir(),function(x)dir(x))
file.list <- unlist(sapply(1:length(files), function(x)return(paste(names(files)[x],files[[x]],sep="/"))))

arrays <- file.list[sapply(cel.list, function(x)grep(x, file.list))]

if(platform==FALSE){

if(length(dataset.list)>1){

dataset.list <- as.factor(dataset.list)

datasets <- tapply(arrays, dataset.list, function(x){

  tmp_affy <- read.affybatch(x)
  tmp.exprs <- mas5(tmp_affy, normalize=TRUE)
  

  
  #return(test)
  return(exprs(tmp.exprs))

})


sizes <- sapply(datasets, function(x)length(unique(rownames(x))))
gs <- tapply(datasets, sizes, function(x)rownames(x[[1]]))

gs.u <- as.factor(unlist(gs))
as.character(sort(unique(gs.u))[tabulate(gs.u)==length(gs)]) -> sel
sel[-grep("AFFX",sel)] -> sel

datasets.sel <- sapply(datasets, function(x)return(x[sel,]))

big_data <- datasets.sel[[1]]
sapply(2:length(datasets.sel), function(x)big_data<<-cbind(big_data, datasets.sel[[x]][,]));



                        }
                  }
else  {

#data.raw.1 <- read.affybatch(arrays)
#data.exprs.1 <- exprs(rma(data.raw.1, normalize=norm))

#dres <- david_assign2(conversion, data.exprs.1)
#dres <- cbind(Genes = rownames(dres), Probes=as.vector(unlist(attributes(dres)[3])), dres)
#dres <- dres[-grep("AFFX",dres[,2]),]

#return(dres)
print('No datasets selected')
}

}

joinSets <- function(set1, set2, other=FALSE){

gn1 <- set1[,1]
gn2 <- set2[,1]

probe1 <- set1[,2]
probe2 <- set2[,2]

sect <- intersect(gn1, gn2)

sel.1 <- match(sect, gn1)
sel.2 <- match(sect, gn2)

if(other==TRUE){
 read.delim("Y:/Normal Muscle Model/Raw files/GSE14901_SM/GSE14901_data_only.txt",sep="\t") -> mat
 if(length(gn1)>length(gn2)){
  mat.final <- mat[sel.1,-1]
 }
 else {
  mat.final <- mat[sel.2,-1]
 }
mat.final <- log2(mat.final)

sam.sel <- as.character(read.delim("Y:/Normal Muscle Model/Raw files/GSE14901_SM/GSE14901_NMM_samples.txt",F)[,1])
mat.final <- mat.final[,match(sam.sel, colnames(mat.final))]

data.final <- cbind(set1[sel.1,], set2[sel.2,-c(1:2)], mat.final)

colnames(data.final)[1] <- "Genes"
colnames(data.final)[2] <- "Probes"

return(data.final)

}

data.final <- cbind(gn1[sel.1,], gn2[sel.2,-c(1:2)])

colnames(data.final)[1] <- "Probes"
colnames(data.final)[2] <- "Genes"

return(data.final)

}
