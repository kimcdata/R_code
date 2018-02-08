################## R FUNCTION TO CONVERT A SET OF TEXT FILES TO A MULTI-SHEET EXCEL FILE ####################

txtToXlsx = function(list_of_files, filename = "merged.files.xlsx", list_of_sheetnames){

require(xlsx)
wb <- createWorkbook(type="xlsx")


for file in list_of_files{

sheet_data = read.delim(file, row.names = 1)
sheet1 <- createSheet(wb, file)

addDataFrame(sheet_data, sheet1)

}

}