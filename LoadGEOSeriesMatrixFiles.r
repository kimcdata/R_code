######################## PARSE GENE EXPRESSION OMNIBUS SERIES MATRIX FILES TO A DATAFRAME ###########################

f = "GSE62564_series_matrix.txt"

LoadGEOSeriesMatrixFile = function(filepath){
  
 f = filepath
 f_lines = readLines(f)
 series_metadata_lines = grep("!Series",f_lines, ignore.case=T)
 print(paste(length(series_metadata_lines), " lines of series metadata",sep=""))
 
 sample_metadata_lines = grep("!Sample",f_lines, ignore.case=T)
 print(paste(length(sample_metadata_lines), " lines of sample metadata", sep=""))
 
 table_start = grep("!series_matrix_table_begin", f_lines, ignore.case=T) + 1
 table_end = grep("!series_matrix_table_end", f_lines, ignore.case=T) - 1
 
 data_lines = table_start:table_end
 print(paste(length(data_lines), " lines of data",sep=""))
 
 #################### SERIES METADATA #######################
 
 series_metadata = f_lines[series_metadata_lines]
 
 series_metadata = strsplit(series_metadata, "\t")
 series_max = max(sapply(series_metadata, length))
 series_metadata_fixed = lapply(series_metadata, function(x)return(c(x, rep(NA, series_max-length(x)))))
 
 series_names = sapply(series_metadata_fixed, function(x)x[1])
 series_names = gsub("!","",series_names)
 
 series_data = sapply(series_metadata_fixed, function(x)x[2])
 series_data = gsub("\"","",series_data)
 
 series_list= as.list(series_data)
 names(series_list) = series_names
 
 #################### SAMPLE METADATA ########################
 
 sample_metadata = f_lines[sample_metadata_lines]
 
 sample_metadata = strsplit(sample_metadata, "\t")
 sample_max = max(sapply(sample_metadata, length))
 sample_metadata_fixed = lapply(sample_metadata, function(x)return(c(x, rep(NA, sample_max-length(x)))))
 
 sample_metadata_names = sapply(sample_metadata_fixed, function(x)x[1])
 sample_metadata_names = gsub("!","",sample_metadata_names)
 
 sample_data = sapply(sample_metadata_fixed, function(x)x[2:length(x)])
 sample_data = gsub("\"","",sample_data)
 
 colnames(sample_data) = sample_metadata_names
 sample_data = data.frame(sample_data, stringsAsFactors = F)
 
 
 ###################### DATA #########################
 
 data = f_lines[data_lines]
  
 data = strsplit(data, "\t")
 
 if(length(data) > 1) {
 
 data_rows = sapply(data, function(x)return(x[1]))
 data_rows = gsub("\"","",data_rows)[-1]

 data_columns = gsub("\"","",data[[1]])[-1]
   
 data_max = max(sapply(data, length))
 
 data_fixed_length = sapply(data, function(x){
   
   return(c(x, rep("",data_max-length(x))))
   
 })
 
 data_fixed_length = t(data_fixed_length[-c(1),-c(1)])
 data_fixed_length = apply(data_fixed_length, 2, as.numeric)
 
 print(paste(ncol(data_fixed_length)," columns of data", sep=""))
 
 rownames(data_fixed_length) = data_rows
 colnames(data_fixed_length) = data_columns
  
 data_fixed_length = data.frame(data_fixed_length)  
 } else {
   
 print("It looks like the file contains no expression data")
 data_fixed_length = "no data to return"
 
 }
 
 ####################### RETURN LIST ####################
 
 final_list = list(series = series_list, sample = sample_data, data = data_fixed_length)
 return(final_list)
 
 
}