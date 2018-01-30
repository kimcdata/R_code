
samfromclip <- function(y){

if(!"package:samr"%in%search()){library(samr)}
x <- read.delim("clipboard",F)
data <- list(x=x, y=y, logged2=TRUE)
samr(data, resp.type="Two class unpaired") -> samr.res
samr.compute.delta.table(samr.res) -> delta.tab
siggenes.table <- samr.compute.siggenes.table(samr.obj=samr.res, data=data, delta.table=delta.tab,all.genes=T)
rbind(siggenes.table$genes.up, siggenes.table$genes.lo) -> res

res[order(as.numeric(res[,"Row"])),] -> res
write.table(res, file=paste("SAM.results.",gsub("\\:","\\.",Sys.time()),".txt",sep=""), sep="\t",quote=FALSE,row.names=FALSE)
}