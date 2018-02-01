gene_expression_file = "~/cbf-fs/CLL_correlation_network/GSE21029_expression_mat.txt"
target_file = "Overlapping.genes.selected.SNPs.txt"
target_pathways_file = "target_pathways.txt"
source("~/R_code/getFunctionalProfile.r")
mat2sif_file = "~/Python/mat2sif.0.2-hubs.py"

corrNetWithFuncAnnot = function(gene_expression_file, target_file, target_pathways_file, python_script_file, fdr = 0.1, results_dir = "results"){

	###################### REQUIRED LIBRARYS #########################

	require(data.table)
	require(clusterProfiler)

	####################### LOAD GENE EXPRESSION MATRIX ###########################

	expr_data = read.delim(gene_expression_file, row.names=1)
	expr_data = expr_data[,grep("Lymph", colnames(expr_data))]

	###################### SUBSET MATRIX FOR TARGET GENES #########################

	target_genes = read.delim(target_file, stringsAsFactors=F)
	target_genes_harm = intersect(rownames(expr_data), target_genes[,1])

	expr_data_subset = expr_data[target_genes_harm,]

	str(expr_data_subset)


	##################### CREATE CORRELATION MATRIX ##############################

	cor_mat = apply(expr_data_subset, 1, function(x){

		result = cor(x, t(expr_data), method="spearman")

	})

	rownames(cor_mat) = rownames(expr_data)

	print("cor mat") 
	str(cor_mat)

	cor_mat_melt = melt(cor_mat, variable.factor = F)

	print("cor mat melt")
	str(cor_mat_melt)

	#################### WRITE CORRELATION MATRIX TO FILE ########################

	write.table(cor_mat, file="correlation_matrix.txt", sep="\t", quote=F, col.names=NA)

	#################### RANDOMISE DATA FOR RESAMPLING PROCEDURE ####################

	cor_rand = cor(t(apply(expr_data, 2, sample)), method="spearman")
	cor_rand[upper.tri(cor_rand, diag = T)] = NA
	cor_rand = melt(cor_rand, na.rm=T)

	print("cor rand")
	str(cor_rand)

	#################### CALCULATE MEAN/SD OF CORRELATION VALUES FROM RANDOMISED DATA #####################

	cor_rand_mean = mean(cor_rand$value)
	cor_rand_sd = sd(cor_rand$value)

	####################### CALCULATE P VALUE MATRIX ##########################

	print("start")

	pos_flag = cor_mat_melt$value > 0
	pos_p = pnorm(cor_mat_melt$value[pos_flag], mean = cor_rand_mean, sd = cor_rand_sd, lower.tail = F)
	neg_p = pnorm(cor_mat_melt$value[!pos_flag], mean = cor_rand_mean , sd = cor_rand_sd, lower.tail = T)

	cor_mat_melt$p_value = rep(1, length(pos_flag))
	cor_mat_melt$p_value[pos_flag] = pos_p
	cor_mat_melt$p_value[!pos_flag] = neg_p



	###################### CORRECT P VALUES FOR MULTIPLE TESTING #######################

	cor_mat_melt$fdr = p.adjust(cor_mat_melt$p_value, "BH")

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

	source("~/R_code/getFunctionalProfile.r")

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
}



