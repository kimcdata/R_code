
corrNetWithFuncAnnot = function(gene_expression_file, target_file, python_script_file, fdr = 0.1){
  
  ####################### LOAD GENE EXPRESSION MATRIX ###########################
  
  expr_data = read.delim(gene_expression_file, row.names=1)

  
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
  
  str(cor_mat)
  
  cor_mat_melt = melt(cor_mat)
  
  str(cor_mat_melt)

  #################### WRITE CORRELATION MATRIX TO FILE ########################
  
  write.table(cor_mat, file="correlation_matrix.txt", sep="\t", quote=F, col.names=NA)
  
  #################### RANDOMISE DATA FOR RESAMPLING PROCEDURE ####################
  
  cor_rand = cor(t(apply(expr_data, 2, sample)), method="spearman")
  cor_rand[upper.tri(cor_rand, diag = T)] = NA
  cor_rand = melt(cor_rand, na.rm=T)
  
  #################### CALCULATE MEAN/SD OF CORRELATION VALUES FROM RANDOMISED DATA #####################
  
  cor_rand_mean = mean(cor_rand$value)
  cor_rand_sd = sd(cor_rand$value)
  
  ####################### CALCULATE P VALUE MATRIX ##########################
  
  pos_flag = cor_mat_melt$value > 0
  pos_p = pnorm(cor_mat_melt$value[pos_flag], mean = cor_rand_mean, sd = cor_rand_sd)
  neg_p = pnrom(cor_mat_melt$value[!pos_flag], mean = cor_rand_mean , sd = cor_rand_sd)
  
  cor_mat_melt$pvalue = rep(1, length(pos_flag))
  cor_mat_melt$p_value[pos_flag] = pos_p
  cor_mat_melt$p_value[!pos_flag] = neg_p
  
  
  
  ###################### CORRECT P VALUES FOR MULTIPLE TESTING #######################
  
  cor_mat_melt$fdr = p.adjust(p_value, "BH")
  
  ###################### STORE ESTIMATE OF CORRELATION VALUES FOR SELECTED FDR VALUES #########################
  
  sig = ifelse(min(cor_mat_melt$fdr) <= fdr, which(cor_mat_melt$fdr < fdr), which.min(cor_mat_melt$fdr))
  corr_thesh = min(abs(cor_mat_melt$value[sig]))
  
  ###################### CREATE NETWORKS FROM CORRELATION MATRIX ##########################
  
  system(paste('python ', python_script_file,' -i correlation_matrix.txt -o correlation_networks.sif.txt -h ', target_file, ' -t ', corr_thresh, sep=""))
  
  ###################### LOAD NETWORKS IN R ###########################
  
  networks = read.delim("correlation_networks.sif.txt", header=F,stringsAsFactors=F)
  
  neighbours = lapply(target_genes_harm, function(x){
    
    r1 = which(networks[,1] == x)
    r2 = which(networks[,2] == x)
    neighbours = unique(c(networks[r1,1],networks[r2,2]))
    return(neighbours)
    
  })
  
  ##################### PERFORM FUNCTIONAL ANALYSIS ON UNIQUE GENES IN EACH NETWORK ########################
  
  functional_profiles = lapply(target_genes_harm, function(x){
    
    getFunctionalProfile(neighbours[[x]], universe = rownames(expr_data), filename = paste("Functional.profile.",x,".txt",sep=""))
    
  })
  
  
  ##################### EXTRACT TARGET PATHWAYS ##########################
  
  target_paathways = read.delim("target_pathways.txt", header=F, stringsAsFactors=F)
  
  sapply(functional_profiles, function(x){
    
    sel = grep(target_pathways, x$Name)
    path = x[sel,]
    
  })
  
}


