
require(data.table)
require(clusterProfiler)
source("~/R_code/getFunctionalProfile.r")

# function to return resampled matrix data in melted format

corrMatResample = function(expr_data, melt_data = FALSE){
	
	#################### RANDOMISE DATA FOR RESAMPLING PROCEDURE ####################

	print("RANDOMISING DATA")
	cor_rand = cor(t(apply(expr_data, 2, sample)), method="spearman")
	
	if(melt_data){
		print("MELTING RANDOMISED CORRELATION MATRIX")
		print("REMOVING UPPER TRIANGLE AND DIAGONAL")
		cor_rand[upper.tri(cor_rand, diag = T)] = NA
		cor_rand = melt(cor_rand, na.rm = T)
		print("MELTED")
		return(cor_rand)
	} else {
		print("CALCULATING MEAN AND SD")
		cor_mean = mean(cor_rand[upper.tri(cor_rand, diag = F)])
		print("COR MEAN")
		print(cor_mean)
		cor_sd = sd(cor_rand[upper.tri(cor_rand, diag = F)])
		print("COR SD")
		print(cor_sd)
		return(c(mean = cor_mean, sd = cor_sd))	
	}
}


# function to return resampled matrix of differential correlation values

diffCorrMatResample = function(expr_data_A, expr_data_B){

#resample data i.e. generate random data
apply(expr_data_A, 2, sample) -> A_resample
apply(expr_data_B, 2, sample) -> B_resample

cor(t(A_resample),method="spearman") -> A_resample_spearman
cor(t(B_resample),method="spearman") -> B_resample_spearman

#create matrix of correlation differences for resampled data
A_resample_spearman - B_resample_spearman -> AB_resample_corrdiff

return(AB_resample_corrdiff)

}

# calculate differential correlation matrix

diffCor = function(data_a, data_b, method="spearman"){

cor_a = cor(t(data_a), method = method)
cor_b = cor(t(data_b), method = method)

diff_cor = cor_a - cor_b

return(diff_cor)

}


# function to return p-values based on mean and sd

corrMatPnorm = function(cor_rand_mean, cor_rand_sd, cor_real_melt, p.adjust.method = "BH"){

	####################### CALCULATE P VALUE MATRIX ##########################

	cor_mat_melt = cor_real_melt
	
	print("start")

	pos_flag = cor_mat_melt$value > 0
	pos_p = pnorm(cor_mat_melt$value[pos_flag], mean = cor_rand_mean, sd = cor_rand_sd, lower.tail = F)
	neg_p = pnorm(cor_mat_melt$value[!pos_flag], mean = cor_rand_mean , sd = cor_rand_sd, lower.tail = T)

	cor_mat_melt$p_value = rep(1, length(pos_flag))
	cor_mat_melt$p_value[pos_flag] = pos_p
	cor_mat_melt$p_value[!pos_flag] = neg_p
	
	cor_mat_melt$fdr = p.adjust(cor_mat_melt$p_value, p.adjust.method)

	return(cor_mat_melt)

}

# fetch a percentile based threshold for correlation

corrPercentileThresh = function(expr_data, perc = 99){

cormat = cor(t(expr_data))
perc_str = paste0(perc,".0%")
thresh = quantile(abs(cormat[upper.tri(cormat, diag=F)]), probs = seq(0,1,length.out=201))[perc_str]
return(thresh)

}


# get node degree from a sif format file

networkDegreeFromSif = function(file, targets){

networks = read.delim(file, header=F,stringsAsFactors=F)

neighbours = sapply(targets, function(x){

	r1 = which(networks[,1] == x)
	r2 = which(networks[,2] == x)
	neighbours = unique(c(networks[r1,2],networks[r2,1]))
	return(neighbours)

})

neighbours_length = sapply(neighbours, length)

return(neighbours_length)

}



corrNetWithFuncAnnot = function(gene_expression_file, target_file, target_pathways_file, python_script_file, fdr = 0.1, alt_thresh = 0.7, results_dir = "results"){

	####################### LOAD GENE EXPRESSION MATRIX ###########################

	expr_data = read.delim(gene_expression_file, row.names=1)

	###################### SUBSET MATRIX FOR TARGET GENES #########################

	target_genes = read.delim(target_file, stringsAsFactors=F)
	target_genes_harm = intersect(rownames(expr_data), target_genes[,1])

	expr_data_subset = expr_data[target_genes_harm,]

	print("READ DATA")
	str(expr_data_subset)

	##################### CREATE CORRELATION MATRIX ##############################

	cor_mat = apply(expr_data_subset, 1, function(x){

		result = cor(x, t(expr_data), method="spearman")

	})

	rownames(cor_mat) = rownames(expr_data)

	print("CORRELATION MATRIX") 
	str(cor_mat)

	cor_mat_melt = melt(cor_mat, variable.factor = F)

	print("CORRELATION MATRIX - MELTED")
	str(cor_mat_melt)

	#################### WRITE CORRELATION MATRIX TO FILE ########################

	write.table(cor_mat, file=paste0(results_dir,"/correlation_matrix.txt"), sep="\t", quote=F, col.names=NA)

	cor_rand_mean_sd = corrMatResample(expr_data, melt_data = FALSE)
		
	#################### CALCULATE MEAN/SD OF CORRELATION VALUES FROM RANDOMISED DATA #####################

	cor_rand_mean = cor_rand_mean_sd["mean"]
	cor_rand_sd = cor_rand_mean_sd["sd"]

	cor_mat_melt = corrMatPnorm(cor_rand_mean, cor_rand_sd, cor_mat_melt)
	

	print("CORRELATION MATRIX WITH STATS - MELTED")
	str(cor_mat_melt)


	###################### STORE ESTIMATE OF CORRELATION VALUES FOR SELECTED FDR VALUES #########################

	sig = ifelse(cor_mat_melt$fdr <= fdr, T, F)

	str(sig)

	minimum_hits = 100

	if(length(which(sig)) < minimum_hits){
		corr_thresh = corrPercentileThresh(expr_data, perc = 99) 
		writeLines(as.character(corr_thresh), con = paste0(results_dir, "/correlation_threshold.txt"))
	} else {
		corr_thresh = min(abs(cor_mat_melt$value[sig]))
	}

	print ("CORRELATION THRESHOLD")

	###################### CREATE NETWORKS FROM CORRELATION MATRIX ##########################

	system(paste0('python ', python_script_file,' -i ', results_dir, '/correlation_matrix.txt -o ', results_dir, '/correlation_networks.sif.txt -t ', corr_thresh))

	###################### LOAD NETWORKS IN R ###########################

	networks = read.delim(paste0(results_dir,"/correlation_networks.sif.txt"), header=F,stringsAsFactors=F)

	print("NETWORKS")
	str(networks)

	neighbours = sapply(target_genes_harm, function(x){

		r1 = which(networks[,1] == x)
		r2 = which(networks[,2] == x)
		neighbours = unique(c(networks[r1,2],networks[r2,1]))
		return(neighbours)

	})

	neighbours_length = sapply(neighbours, length)
	cat("\n\nNEIGHBOURS LENGTH\n\n")
	print(neighbours_length)

	cat("\n\nNEIGHBOURS\n\n")
	str(neighbours)



	##################### PERFORM FUNCTIONAL ANALYSIS ON UNIQUE GENES IN EACH NETWORK ########################

	cat("\nFUNCTIONAL PROFILES STARTED\n\n\n")
	
	functional_profiles = lapply(target_genes_harm, function(x){

		if(length(neighbours[[x]]) > 0){
			cat(paste0("\n\nFUNCTIONAL ANALYSIS GENE: ",x,"\n\n"))
			return(getFunctionalProfile(gene.list = neighbours[[x]], universe = rownames(expr_data), filename = paste0(results_dir,"/Functional.profile.",x,".txt"), organism.db = "org.Hs.eg.db", organism = "human", golevels = F, kegg.organism = "hsa"))
		} else {

			cat(paste0("\n\nGENE: ",x," HAS NO NEIGHBOURS\n\n"))
			res = matrix(NA, ncol = 10, nrow = 1)
			colnames(res) = c("Category","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count")
			return(res)

		}
	})

	save(functional_profiles, file = paste0(results_dir,"/functional.profiles.current.run.RData"))

	cat("\nFUNCTIONAL PROFILES COMPLETED\n\n")
	
	##################### EXTRACT TARGET PATHWAYS ##########################

	target_pathways = as.character(read.delim(target_pathways_file, header=F, stringsAsFactors=F)[,1])

	path_tables = list()
	target_columns = c("pvalue","p.adjust","qvalue","geneID","Count")

	for(p in target_pathways){

		cat(paste0("\nPATHWAY: ",p,"\n\n"))

		path_tables[[p]] = sapply(functional_profiles, function(f){
			str(f)
			search = f$Description == p
			if(any(search)){
				res = f[search,target_columns]
				return(res)

			} else {
				res = rep(NA, 5)
				names(res) = target_columns
				return(res)
			}
		})
	}

	path_tables = lapply(path_tables, function(x)apply(x, 2, unlist))
	names(path_tables) = gsub("[^A-z]",".",names(path_tables))

	for(p in names(path_tables)){

		to_write = t(path_tables[[p]])
		rownames(to_write) = target_genes_harm
		str(to_write)

		write.table(to_write, file = paste0(results_dir,"/Target.pathway.",p,".enrichment.scores.txt"), quote=F, sep = "\t", col.names=NA)

	}
	
	cat("\nPATHWAY TABLES WRITTEN TO FILE\n\n")
	
}


diffCorrNetWithFuncAnnot = function(gene_expression_file_A, gene_expression_file_B, target_file, target_pathways_file, python_script_file, fdr = 0.1, results_dir = "results"){

	####################### LOAD GENE EXPRESSION MATRIX ###########################

	expr_data_A = read.delim(gene_expression_file_A, row.names=1)
	epxr_data_B = read.delim(gene_expression_file_B, row.names=1)

	###################### SUBSET MATRIX FOR TARGET GENES #########################

	target_genes = read.delim(target_file, stringsAsFactors=F)
	target_genes_harm = intersect(rownames(expr_data), target_genes[,1])

	str(expr_data_subset)

	########################### CREATE DIFFERENTIAL CORRELATION MATRIX ##############################
	
	diff_cor_mat = diffCor(expr_data_A, expr_data_B)
	
	cor_mat = diff_cor_mat[target_genes_harm,]

	print("cor mat") 
	str(cor_mat)

	cor_mat_melt = melt(cor_mat, variable.factor = F)

	print("cor mat melt")
	str(cor_mat_melt)

	#################### WRITE CORRELATION MATRIX TO FILE ########################

	write.table(cor_mat, file="correlation_matrix.txt", sep="\t", quote=F, col.names=NA)

	
	#################### CALCULATE MEAN/SD OF CORRELATION VALUES FROM RANDOMISED DATA #####################

	cor_rand = diffCorrMatResample(expr_data_A, expr_data_B, melt = T)
	
	cor_rand_mean = mean(cor_rand$value)
	cor_rand_sd = sd(cor_rand$value)

	cor_mat_melt = corrMatPnorm(cor_rand_mean, cor_rand_sd, cor_mat_melt)

	print("cor mat melt")
	str(cor_mat_melt)


	###################### STORE ESTIMATE OF CORRELATION VALUES FOR SELECTED FDR VALUES #########################

	sig = ifelse(cor_mat_melt$fdr <= fdr, T, F)

	str(sig)

	minimum_hits = 100

	if(length(which(sig)) < minimum_hits){
		corr_thresh = 0.7
	} else {
		corr_thresh = min(abs(cor_mat_melt$value[sig]))
	}

	str(corr_thresh)

	###################### CREATE NETWORKS FROM CORRELATION MATRIX ##########################

	system(paste('python ', python_script_file,' -i correlation_matrix.txt -o correlation_networks.sif.txt -t ', corr_thresh, sep=""))

	###################### LOAD NETWORKS IN R ###########################

	networks = read.delim("correlation_networks.sif.txt", header=F,stringsAsFactors=F)

	print("networks")
	str(networks)

	neighbours = sapply(target_genes_harm, function(x){

		r1 = which(networks[,1] == x)
		r2 = which(networks[,2] == x)
		neighbours = unique(c(networks[r1,2],networks[r2,1]))
		return(neighbours)

	})

	print("neighbours")
	str(neighbours)



	##################### PERFORM FUNCTIONAL ANALYSIS ON UNIQUE GENES IN EACH NETWORK ########################

	functional_profiles = sapply(target_genes_harm, function(x){

		getFunctionalProfile(gene.list = neighbours[[x]], universe = rownames(expr_data), filename = paste0(results_dir,"/Functional.profile.",x,".txt"), organism.db = "org.Hs.eg.db", organism = "human", golevels = F, kegg.organism = "hsa")

	})

	##################### EXTRACT TARGET PATHWAYS ##########################

	target_pathways = as.character(read.delim(target_pathways_file, header=F, stringsAsFactors=F)[,1])

	path_tables = list()
	target_columns = c("pvalue","p.adjust","qvalue","geneID","Count")

	for(p in target_pathways){

		path_tables[[p]] = sapply(functional_profiles, function(f){
			search = f$Description == p
			if(any(search)){
				res = f[search,target_columns]
				return(res)

			} else {
				res = rep(NA, 5)
				names(res) = target_columns
				return(res)
			}
		})
	}

	path_tables = lapply(path_tables, function(x)apply(x, 2, unlist))
	names(path_tables) = gsub("[^A-z]",".",names(path_tables))

	for(p in names(path_tables)){

		to_write = path_tables[[p]]
		str(to_write)

		write.table(t(to_write), file = paste0(results_dir,"/Target.pathway.",p,".enrichment.scores.stxt"), quote=F, sep = "\t", col.names=NA)

	}
	
	return(functional_profiles)
}


