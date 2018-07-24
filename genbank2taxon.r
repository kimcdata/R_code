############# FETCH TAXONOMIC IDS FROM GENBANK IDENTIFIERS (E.G. FROM A BLASTN > NR RESULT) ##############

devtools::install_github("ropensci/taxize")

library('taxize')

########### NCBI API KEY ##############

api_key = '9d2431879162e81f1424da7632cc76202507'

############# LOAD SOME SAMPLE DATA #################

ids = read.delim("C:/dropbox2/CBF/Other/AHollander_MSCs/RNAseq_qc/testblast.ids.out.txt",F,stringsAsFactors = F)

taxize::genbank2uid(id = ids[1], key = api_key)




