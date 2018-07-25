################## FILTER A NETWORK AND COMPUTE JACCARD INDEX BETWEEN NODE NEIGHBOURHOODS #######################

networkJaccard = function(sif, vids, neigh = 3){
 
  require(igraph)

  sif_two_column = sif[,c(1,2)]
  sif_vector = c(t(sif_two_column))
  str(sif_vector)
  
  vids = sort(intersect(sif_vector, vids))
  str(vids)
  
  g = make_graph(sif_vector)
  ns = neighborhood.size(g)
  str(ns)
  ns_sel = which(ns > neigh)
  str(ns_sel)
  ns_vids = vertex.attributes(g)$name[ns_sel]
  str(ns_vids)
  vids = intersect(vids, ns_vids)
  
  s = similarity(g,method = "jaccard", vids = vids)
  rownames(s) = colnames(s) = vids
  return(s)
  
}

thresholdSif = function(sif, thresh, thresh_col, abs = F, mode = "gt"){
 
  mode = tolower(mode)
  if(mode == "gt"){
    if(abs){
      sel = which(abs(sif[,thresh_col]) >= thresh)
      } else if(!abs){
      sel = which(sif[,thresh_col] >= thresh) 
      }
  }
    else if (mode == "lt"){
      if(abs){
        sel = which(abs(sif[,thresh_col]) <= thresh)
      } else if(!abs){
        sel = which(sif[,thresh_col] <= thresh)
      }
    }
  return(sif[sel,])
        
}

networkJaccardHeatmap = function(jaccMat, filter = T, filterValue = 0.5){
 
  require(gplots)
  plot.new()
  
  diag(jaccMat) = 0
  sel = which(apply(jaccMat, 1, max, na.rm = T) > filterValue)
  str(sel)
  jaccMat = jaccMat[sel,sel]
  
  heatmap.2(jaccMat, trace="none", col = bluered(50), mar=c(12,12), hclustfun = function(x)hclust(x, "average"), distfun = function(x)as.dist(1-cor(t(x))))
   
}

############# EXAMPLE ###################

sif = read.delim("Lymph.KEGG.GSEA.network.sif.txt",stringsAsFactors = F)
sif = thresholdSif(sif, 0.01, 3, abs = F, mode = "lt")
sif = thresholdSif(sif, 2, 4, abs = T, mode = "gt")

vids = read.delim("jaccard.index.targets.txt", header = F, stringsAsFactors = F)[,1]
s = networkJaccard(sif, vids)
write.table(s, file = "lymph.network.neighbourhood.jaccard.matrix.txt", sep="\t", quote = F, col.names = NA)
networkJaccardHeatmap(s)
