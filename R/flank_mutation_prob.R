#' predict mutation probility

#' The function to predict mutation probility for each mutation

#' @importFrom data.table fread
#' @param lr_result  trainning result used to predict mutation probability for each genomic position
#' @param inner_pos_anno  the file contains mutations and annotation values surrounded the recurrent position ±5 base pairs
#' @param outer_pos_anno: the file contains mutations and annotation values surrounded the recurrent position ±500 base pairs
#' @param inner_pos_prob: output inner position result
#' @param outer_pos_prob: output outer position result

#' @return None
#' @export







flank_mutation_prob <- function(lr_result, inner_pos_anno, outer_pos_anno, inner_pos_prob, outer_pos_prob )
{
	#Args:
	# - lr_result: trainning result used to predict mutation probability for each genomic position
	# - inner_pos_anno: the file contains mutations and annotation values surrounded the recurrent position ±5 base pairs
	# - outer_pos_anno: the file contains mutations and annotation values surrounded the recurrent position ±500 base pairs
	# - inner_pos_prob: output inner position result
	# - outer_pos_prob: output outer position result
	
	load(lr_result)
	
	inner_chr_info <- fread(inner_pos_anno, showProgress = F, data.table = F)
	names(inner_chr_info) <- c("chr","start","end","bp","reptime","tfbs","gerp","dnase","cpg","gencode")
	inner_chr_info[inner_chr_info == "."] <- 0
	inner_chr_info$reptime <- round(as.numeric(inner_chr_info$reptime), 2)
	inner_chr_info$tfbs <- round(as.numeric(inner_chr_info$tfbs), 2)
	inner_chr_info$gerp <- round(as.numeric(inner_chr_info$gerp), 2)
	inner_chr_info$dnase <- round(as.numeric(inner_chr_info$dnase), 2)
	inner_chr_info$cpg <- round(as.numeric(inner_chr_info$cpg), 2)
	inner_chr_info <- unique(inner_chr_info)

	outer_chr_info <- fread(outer_pos_anno, showProgress = F, data.table = F)
	names(outer_chr_info) <- c("chr","start","end","bp","reptime","tfbs","gerp","dnase","cpg","gencode")
	outer_chr_info[outer_chr_info == "."] <- 0
	outer_chr_info$reptime <- round(as.numeric(outer_chr_info$reptime), 2)
	outer_chr_info$tfbs <- round(as.numeric(outer_chr_info$tfbs), 2)
	outer_chr_info$gerp <- round(as.numeric(outer_chr_info$gerp), 2)
	outer_chr_info$dnase <- round(as.numeric(outer_chr_info$dnase), 2)
	outer_chr_info$cpg <- round(as.numeric(outer_chr_info$cpg), 2)	
	outer_chr_info <- unique(outer_chr_info)

	inner_chr_info$prob <- predict(logi_result, newdata = inner_chr_info, type = "response")
	outer_chr_info$prob <- predict(logi_result, newdata = outer_chr_info, type = "response")	
	
	write.table(inner_chr_info, inner_pos_prob, sep = "\t", row.names = F, quote = F)
	write.table(outer_chr_info, outer_pos_prob, sep = "\t", row.names = F, quote = F)	
}
