##################### Function to load a numeric matrix with non-numeric header rows ###############

# skip == number of header rows before the start of the numeric matrix
# rownames == TRUE/FALSE, does the file contain a column of row descriptors e.g. gene/protein IDs?

loadMixedMatrix = function(file, skip, rownames = T){
  
  if(!rownames){ 
    data = read.delim(file, skip = skip, header = F) 
    header = read.delim(file, nrows = 3, header = F)
  } else {
    data = read.delim(file, skip = skip, row.names = 1, header = F) 
    header = read.delim(file, nrows = 3, row.names = 1, header = F)
  }
  
  tmat= data.frame(t(header),t(data),stringsAsFactors = F)
  return(tmat)
  
}
