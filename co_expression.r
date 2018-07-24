#################### CO-EXPRESSION ANALYSIS ####################

# spearman_matrix = function(df){
# 
# 	res = cor(df, method = "spearman")
# 	return(res)
# 
# }
# 
# pearson_matrix = function(df){
# 
# 	res = cor(df, method = "spearman")
# 	return(res)
# 
# }

# filter the input datasets

gexFilterPerc = function(data, var_perc = 0.2, expr_perc = 0.2){
  
  data_var = apply(data, 1, var)
  data_mean = rowMeans(data)
  
  var_cut = quantile(data_var, probs = var_perc)
  mean_cut = quantile(data_mean, probs = expr_perc)
  
  var_which = which(data_var <= var_cut)
  mean_which = which(data_mean <= mean_cut)
  rows_to_cut = union(var_which, mean_which)
  print(length(var_which))
  print(length(mean_which))
  print(length(rows_to_cut))
  
  res = data[-c(rows_to_cut),]
  return(res)
  
}


# randomise a matrix and set the upper triangle of the co-expression matrix to NA

matrixRandomise = function(df, method = "spearman", melt = T){
  
  cor_rand = cor(t(apply(df, 2, sample)), method=method)
  cor_rand[upper.tri(cor_rand, diag = T)] = NA
  if(melt){
    cor_rand = melt(cor_rand, na.rm = T)
  }
  return(cor_rand)
  
}

# get mean and SD of a randomised co-expression matrix

permutedDist = function(cor_rand){
  
  if(nrow(cor_rand) == ncol(cor_rand)){
    cor_mean = mean(cor_rand, na.rm = T)
    cor_sd = sd(cor_rand, na.rm = T)
  } else {
    cor_mean = mean(cor_rand$value, na.rm = T)
    cor_sd = sd(cor_rand$value, na.rm = T)
    
  }
  return(c(mean = cor_mean, sd = cor_sd))	

}



matrixCoexpressPvalues = function(sif_melt, dist){
  
 m = dist['mean']
 s = dist['sd']

 print(m)
 print(s)
 
 p_mat = as.numeric(rep(NA, length(sif_melt$value)))
 
 str(p_mat)
 
 p_mat[sif_melt$value > 0] = pnorm(sif_melt$value[sif_melt$value > 0], mean = m, sd = s, lower.tail = F)
 p_mat[sif_melt$value < 0] = pnorm(sif_melt$value[sif_melt$value < 0], mean = m, sd = s, lower.tail = T)
 
 return(p_mat)
  
}

##### COUNT EDGES FROM SIF FILE FOR A GIVEN THRESHOLD ####

networkDegreeFromSif = function(file, targets, header = T, threshold = 0, target_column = 3){
  
  networks = read.delim(file, header=header,stringsAsFactors=F)
  
  networks = networks[abs(networks[,target_column]) > threshold,]
  
  neighbours = sapply(targets, function(x){
    
    r1 = which(networks[,1] == x)
    r2 = which(networks[,2] == x)
    
    neighbours = unique(c(networks[r1,2],networks[r2,1]))
    return(neighbours)
    
  })
  
  neighbours_length = sapply(neighbours, length)
  
  return(neighbours_length)
  
}

##### RETURN NEIGHBOURHOOD FOR A GIVEN HUB #####

hubTargetsFromSif = function(file, hub, header = T, threshold = 0.6, target_column = 3){
  
  networks = read.delim(file, header=header,stringsAsFactors=F)
  
  networks = networks[abs(networks[,target_column]) > threshold,]
        
  r1 = which(networks[,1] == hub)
  r2 = which(networks[,2] == hub)
    
  neighbours = unique(c(networks[r1,2],networks[r2,1]))
  return(neighbours)
  
}


coExpressionAnalysisHubs = function(dataset, hubs, analysis_name = "co-expression", outfile, method = "spearman", plotsig = F){
  
  #### CALCULATE CO-EXPRESSION MATRIX BETWEEN HUBS AND ALL OTHER VARIABLES ####
  
  d = dataset
  hub_matches = intersect(hubs, rownames(d))
  d_hub = d[hub_matches,]
  cor_mat = cor(t(d), t(d_hub), method = method)
  
  #### MELT CO-EXPRESSION MATRIX AND WRITE TO SIF FILE ####
  
  m = melt(cor_mat)
  
  #### GENERATE MEAN AND SD FOR RESAMPLED DATA ####
  
  gc()
  perm = permutedDist(cor_rand = matrixRandomise(df = d, method = "spearman"))
  gc()
  
  
  #### GENERATE P-VALUES AND FDR FOR CO-EXPRESSION MATRIX ####
  
  p_val = matrixCoexpressPvalues(sif_melt = m, dist = perm)
  m$p_val = p_val
  m$fdr = p.adjust(m$p_val, "BH")
  
  #### PLOT ABSOLUTE SPEARMAN VERSUS FALSE DISCOVERY RATE
  
  if(plotsig){
    tiff(filename = paste0(outfile,".tiff"),width = 600, height = 600, units = "px")
    plot(abs(m$value), m$fdr, main = analysis_name, ylim=c(0,0.5))
    abline(h = 0.05, lty="dashed", col = "black", lwd = 2)
    abline(h = 0.01, lty="dashed", col = "red", lwd = 2)
    dev.off()
  }
  #### WRITE NETWORK FILE ####
  
  head(m)
  str(m)
  
  write.table(m, file = outfile, sep="\t", quote = F, row.names=F, col.names = T)
  
  return(0)
  
}





############################# OLDER FUNCTIONS FROM DIFFCOR ANALYSIS ###################################

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

##### TEST DATA #####

# data = matrix(rnorm(1000),nrow = 10, ncol = 100)
# 
# cor_mat = cor(t(data))
# p_mat = matrix(999, nrow = nrow(cor_mat), ncol = ncol(cor_mat))
# p_mat[cor_mat > 0] = pnorm(cor_mat[cor_mat > 0], mean = 0, sd = 0.05, lower.tail = F)
# p_mat[cor_mat < 0] = pnorm(cor_mat[cor_mat < 0], mean = 0, sd = 0.05, lower.tail = T)


