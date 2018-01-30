##################### CONVERT REFSEQ IDS TO GENE SYMBOLS #######################

ensembl2symbol = function(ids, bioCpackage = "org.Hs.eg.db"){
  
  require(bioCpackage, character.only = T)
  
  result = select(x = eval(parse(text = bioCpackage)), keys = ids, columns = "SYMBOL", keytype = "REFSEQ")
  
  return(result)
  
}