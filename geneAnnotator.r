######################## Annotate Gene Symbols with Genetic and Pathway Information ######################

geneAnnotator <- function(geneList, chrom = T, cytoband = T, kegg = F, keggfile = "", bp = T, cc = T, mf = T, unique = T, organism = "Human"){
 
  
  organism = tolower(organism)
  
  pkg <- switch(organism, 
         human = c(org = "org.Hs.eg.db", tx = "TxDb.Hsapiens.UCSC.hg38.knownGene"),
         mouse = c(org = "org.Mm.eg.db", tx = "TxDb.Mmusculus.UCSC.mm10.knownGene")
  )
  
  orgCode <- switch(organism, 
                human = "hsa",
                mouse = "mmu"
  )
  
  org <- pkg["org"]
  tx <- pkg["tx"]
  
  sapply(pkg, function(x)require(x, character.only = T))
  require(KEGGREST)
  require(GO.db)
  
  # entrez
  
  entrez <- select(x = eval(parse(text=org)), keys = geneList, columns = "ENTREZID",keytype = "SYMBOL")
  entrez = removeNA(entrez, 2)
  entrezKeys <- entrez$ENTREZID
  print("entrez")
  str(entrez)
  
  write.table(entrez, file = "ANN.entrez.txt", sep="\t", quote=F, col.names=NA)

  # description
  
  desc <- select(x = eval(parse(text=org)), keys = geneList, columns = c("ENTREZID","GENENAME"),keytype = "SYMBOL")
  desc = removeNA(desc, 2)
  print("description")
  str(desc)
  
  write.table(desc, file = "ANN.desc.txt", sep="\t", quote=F, col.names=NA)
  
  # chr
  
  chr <- select(x = eval(parse(text=tx)), keys = entrezKeys, columns = "CDSCHROM",keytype = "GENEID")
  chrUnique <- tapply(chr$CDSCHROM, chr$GENEID, function(x)paste(x, collapse="; "))
  chrUnique <- data.frame(ENTREZID = names(chrUnique), CHR = chrUnique, stringsAsFactors = F)
  print("CHR")
  str(chrUnique)
  
  write.table(chrUnique, file = "ANN.chrunique.txt", sep = "\t", quote=F, col.names=NA)
  
  # cytoband
  if(cytoband){
  cyto <- lapply(entrezKeys, function(x){
    
    return(unlist(mget(as.character(x), org.Hs.egMAP, ifnotfound = NA)))

  })
  
  #print("cytoband")
  #str(cyto)
  
  cyto = sapply(cyto, function(x)paste(x, collapse=", "))
  cyto <- data.frame(ENTREZID = entrezKeys, Cytoband = cyto, stringsAsFactors = F)
  cyto = removeNA(cyto, 1)
    
  print("cytoband")
  str(cyto)
  
  write.table(cyto, file = "ANN.cyto.txt", sep="\t", quote=F, col.names=NA)
  
  } else {
    cyto = "Not run"
  }
  # kegg

  if(kegg){

    kegg_table = read.delim(keggfile, stringsAsFactors = F, header=T)
    kegg_by = by(kegg_table, kegg_table$ENTREZID, function(x)c(unique(x[,1]),paste(x[,2],collapse="; "), paste(x[,3], collapse = "; ")))
    kegg_u = t(sapply(kegg_by, function(x)x))
    colnames(kegg_u) = colnames(kegg_table)
    write.table(kegg_u, file = "ANN.kegg.txt", sep="\t", quote=F, col.names=NA)
    
  } else {
    kegg_u = "Not run"
  }
  
  # bp
  
  go <- select(x = eval(parse(text=org)), keys = geneList, columns = c("ENTREZID","GO"),keytype = "SYMBOL")
  gobp <- go[go$ONTOLOGY == "BP",]
  
  bpAnnot <- select(x = GO.db, keys = gobp$GO, columns = "TERM",keytype = "GOID")
  gobp <- data.frame(gobp, TERM = bpAnnot$TERM,stringsAsFactors = F)
  bpCollapsed <- tapply(gobp$TERM, gobp$ENTREZID, function(x)paste(x, collapse = "; "))
  bpCollapsed <- data.frame(ENTREZID = names(bpCollapsed), GOBP = bpCollapsed, stringsAsFactors = F)
  
  print("GOBP")
  str(bpCollapsed)
  write.table(bpCollapsed, file = "ANN.gobp.txt", sep="\t", quote=F, col.names=NA)
  
  # cc
  
  gocc <- go[go$ONTOLOGY == "CC",]
  ccAnnot <- select(x = GO.db, keys = gocc$GO, columns = "TERM",keytype = "GOID")
  gocc <- data.frame(gocc, TERM = ccAnnot$TERM,stringsAsFactors = F)
  ccCollapsed <- tapply(gocc$TERM, gocc$ENTREZID, function(x)paste(x, collapse = "; "))
  ccCollapsed <- data.frame(ENTREZID = names(ccCollapsed), GOCC = ccCollapsed, stringsAsFactors = F)
  
  print("GOCC")
  str(ccCollapsed)
  write.table(ccCollapsed, file = "ANN.gocc.txt", sep="\t", quote=F, col.names=NA)
  
  # mf
  
  gomf <- go[go$ONTOLOGY == "MF",]
  mfAnnot <- select(x = GO.db, keys = gomf$GO, columns = "TERM",keytype = "GOID")
  gomf <- data.frame(gomf, TERM = mfAnnot$TERM,stringsAsFactors = F)
  mfCollapsed <- tapply(gomf$TERM, gomf$ENTREZID, function(x)paste(x, collapse = "; "))
  mfCollapsed <- data.frame(ENTREZID = names(mfCollapsed), GOMF = mfCollapsed, stringsAsFactors = F)
  
  print("GOMF")
  str(mfCollapsed)
  write.table(mfCollapsed, file = "ANN.gomf.txt", sep="\t", quote=F, col.names=NA)
  
  defaultTable <- merge(x = entrez, y = desc[,-c(2)], by = "SYMBOL")
  print("default table")
  str(defaultTable)
  write.table(defaultTable, file = "ANN.default.table.txt", sep="\t", quote=F, col.names=NA)
  

  listToMerge <- list(defaultTable = defaultTable)

  if(chrom) { listToMerge[["CHR"]] = chrUnique }
  if(cytoband) { listToMerge[["Cytoband"]] = cyto }
  if(kegg) { listToMerge[["kegg"]] = kegg_u }
  if(bp) { listToMerge[["GOBP"]] = bpCollapsed }
  if(cc) { listToMerge[["GOCC"]] = ccCollapsed }
  if(mf) { listToMerge[["GOMF"]] = mfCollapsed }


  print("listToMerge")
  str(listToMerge)

  finalTable <- Reduce(function(x, y) merge(x = x, y = y, all.x = T, by="ENTREZID"), listToMerge)

  #ord <- match(geneList, finalTable$SYMBOL)
  #finalTable[ord,] -> finalTable

  return(finalTable)
  
}


removeNA = function(df, column){
 
  sel = !is.na(df[,column])
  return(df[sel,])
   
}


