#' combine mutations of inner and outer

#' The function to combine the inner mutations and outer mutations for each recurrent mutation

#' @importFrom data.table fread
#' @importFrom reshape cast
#' @param lr_result  list format for the combined result
#' @param inner_pos_anno  the file contains mutations and annotation values surrounded the recurrent position ±5 base pairs
#' @param outer_pos_anno: the file contains mutations and annotation values surrounded the recurrent position ±500 base pairs
#' @param inner_pos_pro: the file contains mutations surrounded the recurrent position ±5 base pairs and their mutation probilities
#' @param outer_pos_pro: the file contains mutations surrounded the recurrent position ±500 base pairs and their mutation probilities

#' @return None
#' @export


clustered_mut_prob <- function(inner_pos_pro, outer_pos_pro, inner_mut, outer_mut, list_result)
{
	#Args:
	# - inner_pos_pro: the file contains mutations surrounded the recurrent position ±5 base pairs and their mutation probilities
	# - outer_pos_pro: the file contains mutations surrounded the recurrent position ±500 base pairs and their mutation probilities
	# - inner_mut: the file contains mutations surrounded the recurrent position ±5 base pairs
	# - outer_mut: the file contains mutations surrounded the recurrent position ±500 base pairs
	# - list_result: list format for the combined result

	inner_prob <- fread(inner_pos_pro, showProgress = F, data.table = F)
	outer_prob <- fread(outer_pos_pro, showProgress = F, data.table = F)
	inner_index <- paste(inner_prob$chr, inner_prob$end, sep = ":")
	outer_index <- paste(outer_prob$chr, outer_prob$end, sep = ":")
	inner_prob_index <- data.frame(inner_prob, "index" = inner_index)
	outer_prob_index <- data.frame(outer_prob, "index" = outer_index)
	
	
	mutation_inner <- read.table(inner_mut,header = F)
	mutation_outer <- read.table(outer_mut, header = F)
	index_inner <- paste(mutation_inner$V8,mutation_inner$V9,mutation_inner$V10,sep = ":")
	index_outer <- paste(mutation_outer$V8,mutation_outer$V9,mutation_outer$V10,sep = ":")
	mutation_inner_index <- data.frame(mutation_inner, index_inner)
	mutation_outer_index <- data.frame(mutation_outer, index_outer)
	unique_index_outer <- unique(index_outer)
	unique_index_inner <- unique(index_inner)
	
	mutation_list_outer <- vector("list", length = length(unique_index_outer))
	mutation_list_inner <- vector("list", length = length(unique_index_inner))
	names(mutation_list_outer) <- unique_index_outer
	names(mutation_list_inner) <- unique_index_inner
	
	for(i in unique_index_inner) mutation_list_inner[[i]] <- mutation_inner_index[mutation_inner_index$index_inner == i,][,1:3]
	for(i in unique_index_outer) mutation_list_outer[[i]] <- mutation_outer_index[mutation_outer_index$index_outer == i,][,1:3]

	
	for(i in unique_index_inner)
	{
		df <- mutation_list_inner[[i]]
		df$count <- 1
		df <- reshape::cast(melt(df, id = c("V1","V2","V3")), V1+V2+V3~variable, sum)
		df_index <- paste(df$V1, df$V3, sep = ":")
		df <- data.frame(df, "index" = df_index)
		mutation_list_inner[[i]] <- merge(inner_prob_index, df, by.x = "index", by.y = "index")	
	}
	
	for(i in unique_index_outer)
	{
		df <- mutation_list_outer[[i]]
		df$count <- 1
		df <- reshape::cast(melt(df, id = c("V1","V2","V3")), V1+V2+V3~variable, sum)
		df_index <- paste(df$V1, df$V3, sep = ":")
		df <- data.frame(df, "index" = df_index)
		mutation_list_outer[[i]] <- merge(outer_prob_index, df, by.x = "index", by.y = "index")	
	}
	
	count_index <- paste(unique_index_inner, unique_index_outer, sep = "-")
	vector_count_index <- vector("list", length = length(unique_index_inner))
	names(vector_count_index) <- count_index
	
	for(i in count_index)
	{
		split1 <- strsplit(i,"-")
		inner <- split1[[1]][1]
		outer <- split1[[1]][2]
		df_inner <- unique(mutation_list_inner[[inner]][,c(2:4,12,16)])
		df_outer <- unique(mutation_list_outer[[outer]][,c(2:4,12,16)])
		vector_count_index[[i]] <- list(df_inner, df_outer)	
	}
	save(vector_count_index, file = list_result)

}








