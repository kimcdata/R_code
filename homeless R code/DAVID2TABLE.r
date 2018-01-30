
DAVID2TABLE <- function(terms){

readClipboard() -> tmp
tmp[-grep("SP_PIR",tmp)] -> tmp
tmp[-grep("UP_SEQ",tmp)] -> tmp
tmp[-grep("Annotation Cluster",tmp)] -> tmp

sapply(sapply(terms, function(x)grep(paste("\t",gsub("\\(","\\\\\\(",x),"\t",sep=""), tmp,ignore.case=T)), function(x)x[1]) -> s
tmp[s] -> tmp
write.table(t(sapply(strsplit(tmp, "\t"), function(x)c(sapply(strsplit(x[2],"_"),function(x)x[2]),x[3],x[6],x[8]))), "clipboard",sep="\t",quote=FALSE,row.names=FALSE,col.names=c("Term","Category","Count","FDR"))

}

