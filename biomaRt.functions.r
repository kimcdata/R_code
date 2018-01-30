####################### BiomaRt functions ########################

library(biomaRt)

ensembl=useMart("ensembl")
#ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)


symbol2x <- function(gl, attributes = c("entrezgene", "hgnc_symbol")){

ensembl = useMart("ensembl",dataset="musmusculus_gene_ensembl")
h2m <- getBM(values = gl, attributes=attributes,mart=ensembl)
return(h2m)

}

affy2symbol <- function(gl, attributes = c("affy_hg_u133_plus_2","hgnc_symbol")){

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
h2m <- getBM(values = gl, attributes=attributes,mart=ensembl)
return(h2m)

}

agilent2symbol <- function(gl, attributes = c("efg_agilent_sureprint_g3_ge_8x60k","with_hgnc")){
  
  ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
  h2m <- getBM(values = gl, attributes=attributes,mart=ensembl)
  return(h2m)
  
}


######################## bioconductor packages ############################

org.Hs.egSYMBOL2EG -> x
as.list(x) -> xx
xx[sapply(xx, length)==1] -> xx
do.call(rbind, xx) -> symbol2eg





