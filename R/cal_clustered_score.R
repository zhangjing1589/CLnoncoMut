#' calculate clustered score

#' The function to calculate clustered score in noncoding regions

#' @param list_result  list format for the combined result
#' @return None
#' @export



cal_clustered_score <- function(list_result)
{
	load(list_result)	

	cal_pois <- function(x, lam)  lam^x*exp(-lam)/prod(1:x)
	
	region_index <- names(vector_count_index)
	chr <- vector(length = length(region_index))
	pos <- vector(length = length(region_index))
	pval <- vector(length = length(region_index))
	
	for(i in 1:length(region_index))
	{
		df_inner <- vector_count_index[[region_index[i]]][[1]]
		df_outer <- vector_count_index[[region_index[i]]][[2]]
		df_inner_index <- paste(df_inner$chr, df_inner$start, df_inner$end, sep = ":")
		df_outer_index <- paste(df_outer$chr, df_outer$start, df_outer$end, sep = ":")
		df_outer_minus_inner <- df_outer[!df_outer_index %in% df_inner_index,]
		df_outer_minus_inner <- df_outer_minus_inner[df_outer_minus_inner$count < 5, ]
		
		if(dim(df_outer_minus_inner)[1] < 1) 
		{
			chr[i] <- strsplit(region_index[i],fixed = F ,"[-:]")[[1]][1]
			pos[i] <- as.numeric(strsplit(region_index[i],fixed = F ,"[-:]")[[1]][3]) - 5
			pval[i] <- Inf
			next
		}
		
		length_inner <- dim(df_inner)[1]
		length_outer <- dim(df_outer_minus_inner)[1]
		
		pvalue_inner <- 1
		pvalue_outer <- 1
		for(j in 1:length_inner)  pvalue_inner <- pvalue_inner*cal_pois(df_inner$count[j],df_inner$prob[j])
		for(j in 1:length_outer)  pvalue_outer <- pvalue_outer*cal_pois(df_outer_minus_inner$count[j], df_outer_minus_inner$prob[j])

		chr[i] <- strsplit(region_index[i],fixed = F ,"[-:]")[[1]][1]
		pos[i] <- as.numeric(strsplit(region_index[i],fixed = F ,"[-:]")[[1]][3]) - 5
		pval[i] <- pvalue_inner/exp(log10(pvalue_outer))
	}
	binding_pvalue <- data.frame(chr, pos, pval)
	write.table(binding_pvalue, "cluster_score.tsv", sep = "\t", row.names = F , quote = F)
}