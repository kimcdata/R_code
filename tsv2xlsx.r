## A function to take an arbitrary number of TSV files and produce an excel workbook ##

tsvToXlsx = function(tsv_files, worksheet_names, xlsx_filename = "workbook.xlsx"){
  
  require(xlsx)
  
  if(missing(worksheet_names)){
    worksheet_names = sapply(tsv_files, function(x)substr(x, 1, 29))
  }
  
  worksheet_names = paste0(1:length(worksheet_names),"_",worksheet_names)
  
  wb<-createWorkbook(type="xlsx")
  
  for(i in 1:length(tsv_files)){
  
  w = worksheet_names[i]
  
  print(w)
  
  f = tsv_files[i]
  sheet_list = list()
  
  #assign(paste0("sheet",i), createSheet(wb, sheetName = w))
  sheet_list$w = createSheet(wb, sheetName = w)
  
  df = read.delim(f, header=T, stringsAsFactors = F, row.names=1)
  df = df[order(df$p.adjust, rownames(df), decreasing = F),]
  rnames = rownames(df)
  df = apply(df, 2, as.character)
  rownames(df) = rnames
  addDataFrame(x = df, sheet = sheet_list$w, col.names = T, row.names = T, showNA = T, characterNA = "NA")
  
  
  }
  
  
  saveWorkbook(wb = wb, file = xlsx_filename)
  
  return("Written to file")
  
}
