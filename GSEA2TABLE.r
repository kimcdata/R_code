

fetchCoreEnrichedGenes <- function(directory, geneset){
 
  fh <- paste(geneset, ".xls",sep="")
  print(fh)
  dirSearch <- dir(directory, pattern=fh)
  print(dirSearch)
  if(length(dirSearch) == 1){
    tab <- read.delim(file.path(directory,fh),header=T,stringsAsFactors = F)
    ceg <- tab$PROBE[tab$CORE.ENRICHMENT == "Yes"]
    ceg <- paste(ceg, collapse=",")
    print("hit")
    return(ceg)
  } else {
    print("swing and a miss")
    return("")
  }  
}

sapply(gseaRuns, function(x){
  sapply(target.regulons, function(y){
    
    fetchCoreEnrichedGenes(paste("aureus/gsea/",x,sep=""),y)
    
  })
}) -> union.core.enriched

sapply(target.regulons, function(x)unique(unlist(union.core.enriched[x,]))) -> regulon.core.enriched
sapply(regulon.core.enriched, function(x)x[x!=""]) -> regulon.core.enriched


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
  
# convert each GSEA result into a matrix
  
lapply(gsea.res, function(x){
  pos.core.enr.vector <- rep("",length(x$nes.pos))
  names(pos.core.enr.vector) <- names(x$nes.pos)
  neg.core.enr.vector <- rep("", length(x$nes.neg))
  names(neg.core.enr.vector) <- names(x$nes.neg)
  pos.core.enr.vector[names(x$pos.core)] <- x$pos.core
  neg.core.enr.vector[names(x$neg.core)] <- x$neg.core
  cbind(c(x$nes.pos, x$nes.neg),c(x$fdrq.pos, x$fdrq.neg),c(pos.core.enr.vector,neg.core.enr.vector))
}) -> gsea.res.mat
  
  
# sort each GSEA result by pathway
lapply(gsea.res.mat, function(x)x[order(rownames(x)),]) -> gsea.res.mat

# create column headings for the final data matrix
rep(c("nes","fdr","enr"),length(tar)) -> ty
paste(ty,names(tar)[expand.grid(1:3,1:length(tar))[,2]]) -> gsea.cols

# combine the gsea results into one matrix
do.call(cbind, gsea.res.mat) -> gsea.res.mat
colnames(gsea.res.mat) <- gsea.cols

# write the results to a file
write.table(gsea.res.mat, file="gsea.merged.txt", sep="\t", quote=FALSE, col.names=NA)
  
#### write a gmx file of the core enriched genes #####

core.list <- lapply(gsea.res, function(x){
  
  l <- list(POS=x$pos.core, NEG=x$neg.core)
  sapply(l, function(y)unlist(strsplit(paste(y, collapse=", "),", ")))
    
  })

max(sapply(core.list, function(x)sapply(x, length))) -> m
lapply(core.list, function(x)sapply(x, function(y)c(y, rep("", m-length(y))))) -> enr.mat
do.call(cbind, enr.mat) -> enr.mat
colnames(enr.mat) <- paste(names(core.list)[sort(rep(1:length(tar),2))], colnames(enr.mat), sep="-")
write.table(enr.mat, file="E.Coli.Union.core.enriched.genes.GSEA.gmx.txt", sep="\t", quote=FALSE, row.names=FALSE)



######################### GSEA NES SCORE HEATMAP #######################

library(gplots)

gsea.res.mat[,grep("nes",colnames(gsea.res.mat))] -> gsea.nes
gsea.res.mat[,grep("fdr",colnames(gsea.res.mat))] -> gsea.fdr
gsub(".FC.*","",colnames(gsea.nes)) -> colnames(gsea.nes)
dimnames(gsea.nes) -> dimn
apply(gsea.nes, 2, as.numeric) -> gsea.nes
dimnames(gsea.nes) <- dimn

dimnames(gsea.fdr) -> dimn
apply(gsea.fdr, 2, as.numeric) -> gsea.fdr
dimnames(gsea.fdr) <- dimn

apply(gsea.nes, 1, function(x)max(abs(x))) -> gsea.nes.max
which(gsea.nes.max > 1.9) -> sel

apply(gsea.fdr, 1, function(x)min(abs(x))) -> gsea.fdr.min
which(gsea.fdr.min < 0.05) -> sel

gsea.nes[sel,] -> gsea.hm
gsea.fdr[sel,] -> gsea.cellnote
apply(gsea.cellnote, 2, as.numeric) -> gsea.cellnote

x11()
heatmap.2(gsea.hm, trace="none", col=bluered(50), mar=c(20,10),cellnote = round(gsea.cellnote,2), notecol = "white")








