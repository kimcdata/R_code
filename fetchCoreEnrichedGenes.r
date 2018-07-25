fetchCoreEnrichedGenes <- function(directory, geneset){
 
  fh <- paste(geneset, ".xls",sep="")
  #print(fh)
  dirSearch <- dir(directory, pattern=fh)
  #print(dirSearch)
  if(fh != dirSearch) { print("MISMATCH BETWEEN FILENAME AND FILE FOUND") }
  if(length(dirSearch) == 1){
    tab <- read.delim(file.path(directory,fh),header=T,stringsAsFactors = F)
    ceg <- tab$PROBE[tab$CORE.ENRICHMENT == "Yes"]
    ceg <- paste(ceg, collapse=",")
    #print("hit")
    return(ceg)
  } else {
    print("swing and a miss")
    return("")
  }  
}