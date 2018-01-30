# bioconductor annotation code #

# annotate a list of genes with entrez IDs and KEGG pathways

library(KEGG.db)
library(org.Mm.eg.db)


columns(org.Mm.eg.db)

gl <- readClipboard()

select(org.Mm.eg.db, gl, "PATH",keytype="SYMBOL") -> path.sel

KEGGPATHID2NAME -> x
mappedkeys(x) -> keys
as.list(x[keys]) -> keggpath2name
cbind(KEGGID = names(keggpath2name),KEGGNAME = unlist(keggpath2name)) -> keggpath2name

org.Mm.egPATH -> x
mappedkeys(x) -> sel
as.list(x[sel]) -> entrez2kegg

cbind(ENTREZ = rep(names(entrez2kegg),sapply(entrez2kegg, length)), KEGGID = unlist(entrez2kegg)) -> entrez2kegg.mat

org.Mm.egALIAS2EG -> x
mappedkeys(x) -> sel
as.list(x[sel]) -> sym2eg

cbind(SYMBOL = rep(names(sym2eg),sapply(sym2eg, length)), ENTREZ = unlist(sym2eg)) -> sym2eg.mat

merge(x = entrez2kegg.mat, y = keggpath2name, all.x = TRUE, by.x = "KEGGID", by.y = "KEGGID", stringsAsFactors=FALSE,drop=FALSE) -> t1

apply(t1, 2, as.character) -> t1

tapply(t1[,"KEGGNAME"], t1[,"ENTREZ"], function(x)x) -> t1.red

cbind(ENTREZ =names(t1.red), KEGGNAME = sapply(t1.red, function(x)paste(x, collapse="; "))) -> t1.red

merge(x = t1.red, y = sym2eg.mat, by.x = "ENTREZ", by.y="ENTREZ") -> t2
apply(t2, 2, as.character) -> t2

tapply(t2[,"KEGGNAME"], t2[,"SYMBOL"], function(x)x) -> sym2kegg

#by(t2, t2[,"SYMBOL"], function(x)c(unique(x[,"SYMBOL"]),paste(x[,"ENTREZ"],collapse=" ;"),paste(x[,"KEGGNAME"],collapse=" ;"))) -> t2.by

match(gl, t2[,"SYMBOL"]) -> sel

writeClipboard(t2[sel,"KEGGNAME"])

#######################################################################