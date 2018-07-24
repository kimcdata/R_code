
combineGSEARuns <- function(directory, thresh=0.05, nes = T, fdr = TRUE, coreEnrichedGenes = T){
  
  getwd() -> starting_dir
  on.exit(setwd(starting_dir))
  setwd(directory)
  
  dir(pattern="*GseaPreranked", full.names = T) -> tar
  n_runs <- length(tar)
  names(tar) <- unlist(sapply(strsplit(tar, "/"), function(x)x[2]))

  lapply(tar, function(i){

    curr_run <- i
    dir(curr_run) -> file.inv
    print(i)

    
    print("reading files")
    
    #file.path(curr_run, "gsea_report_for_na_pos
    
    # read the GSEA result files
    
    neg <- read.delim(file.path(curr_run,grep("gsea_report_for_na_neg.*xls", dir(curr_run), value=T)), stringsAsFactors = F)
    rownames(neg) <- neg$NAME
    pos <- read.delim(file.path(curr_run,grep("gsea_report_for_na_pos.*xls", dir(curr_run), value=T)), stringsAsFactors = F)
    rownames(pos) <- pos$NAME
    
    listNegCoreEnrichedGenes <- c()
    sapply(neg$NAME, function(x){
      listNegCoreEnrichedGenes <<- c(listNegCoreEnrichedGenes,fetchCoreEnrichedGenes(directory = i, geneset = x))
    })
    listPosCoreEnrichedGenes <- c()
    sapply(pos$NAME, function(x){
      listPosCoreEnrichedGenes <<- c(listPosCoreEnrichedGenes,fetchCoreEnrichedGenes(directory = i, geneset = x))
    })
    

    gsea.res <- list(nes ="", fdr="", coreEnrichedGenes="")
    
    if(nes){
      gsea.res$nes <- list(nes.pos = pos$NES, nes.neg = neg$NES)
      names(gsea.res$nes$nes.pos) <- rownames(pos)
      names(gsea.res$nes$nes.neg) <- rownames(neg)
    }
    if(fdr){
      gsea.res$fdr <- list(fdr.pos = pos$FDR.q.val, fdr.neg = neg$FDR.q.val)
      names(gsea.res$fdr$fdr.pos) <- rownames(pos)
      names(gsea.res$fdr$fdr.neg) <- rownames(neg)
    }
    if(coreEnrichedGenes){
      gsea.res$coreEnrichedGenes <- list(core.pos = listPosCoreEnrichedGenes, core.neg = listNegCoreEnrichedGenes)
      names(gsea.res$coreEnrichedGenes$core.pos) <- rownames(pos)
      names(gsea.res$coreEnrichedGenes$core.neg) <- rownames(neg)
    }
    
    return(gsea.res)

  }) -> gseaSummary

  print(names(gseaSummary))
  return(gseaSummary)
  
}