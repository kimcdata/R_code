############## Linux bash commands to download the KEGG files ##################

# fetch xml files for all pathways
curl -s "http://rest.kegg.jp/list/pathway/T01001" | cut -f 1 | while read A; do curl -o "${A}.xml" "http://rest.kegg.jp/get/${A}/kgml" ; done

#rename the files to remove the ":" and replace with "."
rename 's/:/./g' *
  

############## START R CODE #################

### Installing and loading the KEGGgraph library to convert KEGG pathways to gene networks

source("https://bioconductor.org/biocLite.R")
biocLite("KEGGgraph")
library(KEGGgraph)

### Create a new file to contain the KEGG networks

f <- file(description = "KEGG_sif.txt",open = "w")
f_ann = file(description = "KEGG.annotation.txt", open = "w")

### Store a variable with the file paths of the KEGG pathway XML files

files <- dir("keggxml/",full.names = T)

### Start to loop through the files, reading each one and converting it to a network format

sapply(files, function(x){

p <- x

### convert the KEGG pathway to a network format
as.matrix(parseKGML2DataFrame(file = p,reactions = F)) -> kp

### check the pathway contains at least 1 connection (edge)
if(nrow(kp)>0){
  
  # a bit of cleanup to remove unwanted text
  kp[,1] <- gsub("hsa:","",kp[,1])
  kp[,2] <- gsub("hsa:","",kp[,2])
  
  # fetch the KEGG pathway ID from the filename
  pathid <- unlist(strsplit(p, "\\."))[2]
  # add the KEGG pathway ID to the network edges
  cbind(kp, rep(pathid, nrow(kp))) -> kp
  
  # add the KEGG pathway edges to the file
  write.table(kp, f, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
})

### close the file after all the KEGG pathways have been converted to network format and saved
close(f)

### Convert the KEGG pathways to a 3 column annotation file format (entrez id, symbol, pathway)

sapply(files, function(x){
  
  p <- x
  
  ### fetch pathway title and number
  k = parseKGML(file = p)
  title = k@pathwayInfo@title
  pathno = k@pathwayInfo@number
  
  ### convert the KEGG pathway to a network format
  as.matrix(parseKGML2DataFrame(file = p,reactions = F)) -> kp
  
  ### check the pathway contains at least 1 connection (edge)
  if(nrow(kp)>0){
    
    # a bit of cleanup to remove unwanted text
    kp[,1] <- gsub("hsa:","",kp[,1])
    kp[,2] <- gsub("hsa:","",kp[,2])
    
    unique_genes = unique(c(kp[,1],kp[,2]))
    path_table = data.frame(ENTREZID = unique_genes, PATH_NO = rep(pathno, length(unique_genes)), PATH_NAME = rep(title, length(unique_genes)))
    
    # add the KEGG pathway edges to the file
    write.table(path_table, f_ann, sep="\t", quote=FALSE, row.names=FALSE, col.names=T)
  }
})

### close the file after all the KEGG pathways have been converted to network format and saved
close(f_ann)







################## PATHWAY ANNOTATION (NOT GENE-GENE CONNECTIONS, EVERY GENE IN EVERY PATHWAY) ################

f_ann = file(description = "KEGG.annotation.txt", open = "w")

### Convert the KEGG pathways to a 3 column annotation file format (entrez id, symbol, pathway)

sapply(files, function(x){
  
  p <- x
  
  ### fetch pathway title and number
  k = parseKGML(file = p)
  title = k@pathwayInfo@title
  pathno = k@pathwayInfo@number
  
  ### convert the KEGG pathway to a network format
  path_genes = unique(unlist(sapply(k@nodes, function(y)y@name)[sapply(k@nodes, function(x)x@type) == "gene"]))
  path_genes = gsub("hsa:","",path_genes)
  
  ### check the pathway contains at least 1 connection (edge)
  if(length(path_genes)>0){

    path_table = data.frame(ENTREZID = path_genes, PATH_NO = rep(pathno, length(path_genes)), PATH_NAME = rep(title, length(path_genes)))
    
    # add the KEGG pathway edges to the file
    write.table(path_table, f_ann, sep="\t", quote=FALSE, row.names=FALSE, col.names=T)
  }
})

### close the file after all the KEGG pathways have been converted to network format and saved
close(f_ann)



############### There will be duplicate edges in the KEGG network, we should consolidate them ################

library(clusterProfiler)

### load the network file we created

sif <- read.delim("KEGG_sif.txt",F,stringsAsFactors = FALSE,colClasses = "character")

### use the bitr() function within the clusterProfiler library to convert the NCBI gene IDs to symbols
### write a  file containing the ID translations that we will use in Cytoscape

entrez2symbol <- bitr(geneID = c(sif[,1],sif[,2]),fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
write.table(entrez2symbol, file="entrez2symbol.for.KEGG.network.nodes.txt",sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)

### Manipulate the KEGG network (sif) so that it only contains unique gene-gene connections, collapse the relationship types and the pathways

apply(sif, 1, function(x)paste(sort(c(x[1],x[2])),collapse="_")) -> sif_labels
by(sif, sif_labels, function(x){
  lab <- sort(c(x[1,1],x[1,2]))
  c3 <- paste(unique(x[,3]),collapse=", ")
  c4 <- paste(unique(x[,4]),collapse=", ")
  return(c(lab,c3,c4))
  
  }) -> sif.by

t(sapply(sif.by, function(x)x)) -> sif.u

### Add some column headings and write the network to a file

colnames(sif.u) <- c("GeneA","GeneB","InteractionTypes","Pathways")
write.table(sif.u, file="KEGG_sif_unique.txt",sep="\t",quote=FALSE, row.names=FALSE)