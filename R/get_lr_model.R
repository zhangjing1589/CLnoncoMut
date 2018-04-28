#' Training logistic regression model 

#' The function to train logistic regression model of annotated genomic regions

#' @importFrom data.table fread
#' @param pos_annotation  the file contains genomic base pair features: base pair(bp), replicating time(reptime), TFBS signal(tfbs), GERP RS score(gerp), DNaseI signal(dnase), CpG Island(cpg), transcript status(gencode_anno) 
#' @param pos_mut  the file contains somatic mutations(single base substitution) in ICGC database 
#' @param lr_result trainning result used to predict mutation probability for each genomic position
#' @return None
#' @export


get_lr_model <- function(pos_annotation, pos_mut, lr_result)
{	
	#Args:
	# - pos_annotation: the file contains genomic base pair features: base pair(bp), replicating time(reptime), TFBS signal(tfbs), GERP RS score(gerp), DNaseI signal(dnase), CpG Island(cpg), transcript status(gencode_anno) 
	# - pos_mut: the file contains somatic mutations(single base substitution) in ICGC database 
	# - lr_result: trainning result used to predict mutation probability for each genomic position
	
	all_anno <- fread(pos_annotation, showProgress = F, data.table = F)
	names(all_anno) <- c("chr","start","end","bp","reptime","tfbs","gerp","dnase","cpg","gencode")
	all_mutation <- fread(pos_mut, showProgress = F, data.table = F)
	names(all_mutation) <- c("chr","start","end","ref","alt","donor","mutid")

	# match mutations to the annotated regions 
	all_mutation_index <- data.frame(all_mutation[,c(1:3,6)], "count" = 1)
	match_mut <- merge(all_anno, all_mutation_index, by.x = c("chr","start","end"), by.y = c("chr","start","end"), all.x = T)
	match_mut$count[is.na(match_mut$count)] <- 0
	match_mut$bp <- as.factor(match_mut$bp)
	match_mut$gencode <- as.factor(match_mut$gencode)
		
	logi_result <- glm(count ~ bp + reptime + tfbs + gerp + dnase + cpg + gencode, family = binomial(), data = match_mut)
	
	save(logi_result, file = lr_result)
}
