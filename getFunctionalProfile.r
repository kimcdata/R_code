# Gene Ontology and KEGG pathway enrichment in R using the clusterProfiler library
#
# Requires packages 'clusterProfiler', 'GO.db' and the bioconductor annotation library for your species, e.g. 'org.Mm.eg.db' for mouse to be installed
#
# variables:
# gene.list = your target genes for enrichment analysis as Gene Symbols
# universe = the entire list of EXPRESSED genes in your dataset
# organism = one of the following [CASE SENSITIVE!]:
#
# "anopheles" 
# "arabidopsis" 
# "bovine" 
# "canine"
# "chicken" 
# "chimp" 
# "coeli-color" 
# "ecolik12"
# "ecsakai" 
# "fly" 
# "gondii"
# "human" 
# "malaria" 
# "mouse" 
# "pig"
# "rat"
# "rhesus" 
# "worm" 
# "xenopus" 
# "yeast"
# "zebrafish"
#
# organism.db = the bioconductor library containing annotation information for your species, mouse = "org.Mm.eg.db", human = "org.Hs.eg.db", zebrafish = "org.Dr.eg.db", worm (elegans) = "org.Ce.eg.db"
# golevels = TRUE/FALSE, if true the distance (level) from each GO term to the "root" of the Gene Ontology tree will be calculated, will take 5-10 minutes for a standard size gene list (500-1000 genes)
# writeTable = TRUE/FALSE, if true the data will be written to a tab delimitted text file
# filename = the name of the file written when writeTable = TRUE
# pvalue = threshold for selecting significant GO and KEGG terms (adjusted for multiple testing)
# qvalue = False Discovery Rate (FDR) threshold for selecting GO and KEGG terms
# p.adjust.method = method for adjusting p-values for multiple testing, the same options as p.adjust()
# kegg = TRUE/FALSE, if true kegg pathway enrichment will try to be done, but can cause errors with some types of ID
# kegg.organism = organism code recognised by the kegg database e.g. mmu, hsa
#
# To retreive the full Gene Ontology and KEGG analysis, leave both 'pvalue' and 'qvalue' as 1.

#set organism
#organism <- "human"
#organism.db <- "org.Hs.eg.db"

#set organism
#organism <- "mouse"
#organism.db <- "org.Mm.eg.db"

#set organism
#organism <- "ecolik12"
#organism.db <- "org.EcK12.eg.db"

#set filename for output
#filename <- "yourFunctionalProfile.output.txt"

#example function call (golevels = FALSE for fast run):
#  getFunctionalProfile(gene.list = mouse.up, 
#					universe = mouse.universe, 
#					organism=organism, 
#					organism.db=organism.db, 
#					golevels=FALSE,
#					writeTable=TRUE,
#					filename="mouse.up.full.functional.profile.txt")

getFunctionalProfile <- function(gene.list, universe, organism = organism, organism.db = organism.db, golevels=TRUE,writeTable=TRUE, filename=filename, pvalue = 1, qvalue = 1, p.adjust.method="BH",kegg.organism="hsa",id_type = "symbol", minGSSize = 1, maxGSSize = 50000, kegg=T, kegg2symbol=T){
  
  library(clusterProfiler)
  library(GO.db)
  eval(parse(text=paste("require(",organism.db,")",sep="")))

  if(id_type == "symbol"){

    print("convering symbols to entrez IDs")    

    gene = bitr(gene.list, fromType="SYMBOL",toType="ENTREZID",OrgDb=organism.db)[,2]
    str(gene)
  } else if(id_type == "entrez"){
    gene = gene.list
    
  } else {
    print("please use 'symbol' or 'entrez' for id_type")
    return(0)
  }
  
  print(kegg)
  
  if(!missing(universe)){
    if(id_type == "symbol"){

      print("converting universe symbols to entrez IDs")
      universe = bitr(universe, fromType="SYMBOL",toType="ENTREZID",OrgDb=organism.db)[,2]
      str(universe)
    } else if(id_type == "entrez"){
     universe = universe
    } else {
      print("please use 'symbol' or 'entrez' for id_type of gene list and universe")
      return(0)
    }

    print("BP")
    # get Biological Process (BP) GO term enrichment
    enrichGO(gene = gene, 
             universe=universe, 
             OrgDb=organism.db,
             ont="BP",
             pAdjustMethod=p.adjust.method, 
             pvalueCutoff=pvalue, 
             qvalueCutoff=qvalue,
             readable=TRUE,
             minGSSize = minGSSize,
             maxGSSize = maxGSSize) -> ego_bp
    
    # get Cellular Component (CC) GO term enrichment
    print("CC")
    enrichGO(gene = gene, 
             universe=universe, 
             OrgDb=organism.db,
             ont="CC",
             pAdjustMethod=p.adjust.method, 
             pvalueCutoff=pvalue, 
             qvalueCutoff=qvalue,
             readable=TRUE,
             minGSSize = minGSSize,
             maxGSSize = maxGSSize) -> ego_cc
    
    # get Molecular Function (MF) GO term enrichment
    print("MF")
    enrichGO(gene = gene, 
             universe=universe, 
             OrgDb=organism.db,
             ont="MF",
             pAdjustMethod=p.adjust.method, 
             pvalueCutoff=pvalue, 
             qvalueCutoff=qvalue,
             readable=TRUE,
             minGSSize = minGSSize,
             maxGSSize = maxGSSize) -> ego_mf
    
    # get KEGG term enrichment
    if(kegg){
      enrichKEGG(gene = gene,
               universe=universe, 
               organism=kegg.organism,
               pvalueCutoff = pvalue,
               qvalueCutoff = qvalue,
               use_internal_data = FALSE,
               minGSSize = minGSSize,
               maxGSSize = maxGSSize,
               pAdjustMethod = "BH") -> ekegg
    }
    
  } else { 
    
    # get Biological Process (BP) GO term enrichment
    enrichGO(gene = gene, 
             OrgDb=organism.db,
             ont="BP",
             pAdjustMethod=p.adjust.method, 
             pvalueCutoff=pvalue, 
             qvalueCutoff=qvalue,
             readable=TRUE,
             minGSSize = minGSSize,
             maxGSSize = maxGSSize) -> ego_bp
    
    # get Cellular Component (CC) GO term enrichment
    enrichGO(gene = gene, 
             OrgDb=organism.db,
             ont="CC",
             pAdjustMethod=p.adjust.method, 
             pvalueCutoff=pvalue, 
             qvalueCutoff=qvalue,
             readable=TRUE,
             minGSSize = minGSSize,
             maxGSSize = maxGSSize) -> ego_cc
    
    # get Molecular Function (MF) GO term enrichment
    enrichGO(gene = gene, 
             OrgDb=organism.db,
             ont="MF",
             pAdjustMethod=p.adjust.method, 
             pvalueCutoff=pvalue, 
             qvalueCutoff=qvalue,
             readable=TRUE,
             minGSSize = minGSSize,
             maxGSSize = maxGSSize) -> ego_mf
    
    # get KEGG term enrichment
    if(kegg){
    enrichKEGG(gene = gene,
               organism=kegg.organism,
               pvalueCutoff = pvalue,
               qvalueCutoff = qvalue,
               use_internal_data = FALSE,
               minGSSize = minGSSize,
               maxGSSize = maxGSSize) -> ekegg
    }
    
  }
  
  
  ego_bp@result -> bp
  cbind(Category = as.character(rep("BP",nrow(bp))),bp) -> bp
  
  ego_cc@result -> cc
  cbind(Category = as.character(rep("CC",nrow(cc))),cc) -> cc
  
  ego_mf@result -> mf
  cbind(Category = as.character(rep("MF",nrow(mf))),mf) -> mf
  
  if(kegg){
    ekegg@result -> ekegg
    cbind(Category = rep("KEGG",nrow(ekegg)),ekegg) -> ekegg
  
    if(kegg2symbol){
		print("converting KEGG IDs back to symbols")
		ekegg$geneID <- sapply(ekegg$geneID, function(x){
  	    glk <- unlist(strsplit(x, "/"))

        
        gls <- bitr(geneID = glk,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = organism.db,drop = T)[,2]
        gll <- paste(gls, collapse = "/")
      })
      
    }
  }
  
  if(kegg){  
    rbind(bp,cc,mf,ekegg) -> full.table 
  } else {
    rbind(bp,cc,mf) -> full.table 
  }
  
  gr.c <- 4
  bgr.c <- 5
  
  if(golevels){
    
    gr.c <- 5
    bgr.c <- 6
    
    apply(full.table, 1, function(x){
      if(x[1]!="KEGG"){
        return(distance_to_root(goterm=x[2],ont=x[1]))
      } else { return("KEGG") }
    }) -> level
    cbind(GORootDistance=level, full.table) -> full.table
  }
  
  #cat(colnames(full.table),file="test.txt",sep=" ")
  full.table[,gr.c] -> gr
  full.table[,bgr.c] -> bgr
  
  grn <- sapply(gr, function(x)eval(parse(text=x)))
  bgrn <- sapply(bgr, function(x)eval(parse(text=x))) 
  
  full.table[,gr.c] <- grn
  full.table[,bgr.c] <- bgrn
    
  if(writeTable){ write.table(as.matrix(full.table), file=filename,sep="\t",quote=FALSE,row.names=FALSE) }
  
  return(full.table)
  
  
}

distance_to_root <- function(goterm, ont, jumps=0){
  
  switch(ont,
         BP = { parent <- as.character(unlist(mget(goterm, GOBPPARENTS))[1])},
         CC = { parent <- as.character(unlist(mget(goterm, GOCCPARENTS))[1])},
         MF = { parent <- as.character(unlist(mget(goterm, GOMFPARENTS))[1])}
  )
  if(parent=="all"){
    return(jumps)
  } 
  else { 
    jumps <- jumps+1
    distance_to_root(unlist(parent),ont=ont,jumps=jumps)
  }
  
}



bitr2 = function (geneID, fromType, toType, OrgDb, drop = TRUE) 
{
    idTypes <- idType(OrgDb)
    msg <- paste0("should be one of ", paste(idTypes, collapse = ", "), 
        ".")
    if (!fromType %in% idTypes) {
        stop("'fromType' ", msg)
    }
    if (!all(toType %in% idTypes)) {
        stop("'toType' ", msg)
    }
    geneID %<>% as.character %>% unique
    db <- load_OrgDb(OrgDb)
    res <- suppressMessages(select(db, keys = geneID, keytype = fromType, 
        columns = c(fromType, toType)))
    ii <- which(is.na(res[, 2]))
    if (length(ii)) {
        n <- res[ii, 1] %>% unique %>% length
        if (n) {
            warning(paste0(round(n/length(geneID) * 100, 2), 
                "%"), " of input gene IDs are fail to map...")
        }
        if (drop) {
            res <- res[-ii, ]
        }
    }
    return(res)
}
