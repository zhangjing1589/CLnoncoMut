#CLnoncoMut: Clustered Noncoding Mutations

This package is used to prioritize functional noncoding mutations in cancer genomes

## Running CLnoncoMut ##

1. Install the package from github
> library(devtools)  
> install_github("zhangjing1589/CLnoncoMut")

2. Download test file :  
	"all_annotation_chr1.tsv": This file contains the values of base pair(bp), replicating time(reptime), TFBS signal(tfbs), GERP RS score(gerp), DNaseI signal(dnase), CpG Island(cpg), transcript status(gencode_anno) for each base on chromosome 1, any base with NA value was removed.  
	"all_mut_chr1.bed"  This file contains the somatic mutations(single base substitution) on chromosome 1 in ICGC database for 4881 donors, minor allele frequency >= 0.05 were removed.   
	"inner_mut_chr1.bed"; "outer_mut_chr1.bed"  These two files contains mutations closed to each recurrent mutation(>= 5). File 'inner_mut_chr1.bed' contains mutations surrounded the recurrent position ¡À5 base pairs, file 'outer_mut_chr1.bed' contains mutations surrounded the recurrent position ¡À500 base pairs   
	"inner_anno_chr1.tsv"; "outer_anno_chr1.tsv"  These two files contains values of bp, reptime, tfbs, gerp, dnase, cpg, gencode_anno for inner_mut_chr1.bed and outer_mut_chr1.bed  
	These files can be found in R package installing library : /path to R package install/CLnoncoMut/data/  or downloaded from  https://github.com/zhangjing1589/CLnoncoMut/tree/master/data

3. Run get_lr_model() to train logistic regression model of annotated genomic regions
	get_lr_model("all_annotation_chr1.tsv","all_mut_chr1.bed","lr_result.Rdata")
	
4. Run flank_mutation_prob() to predict mutation probility for each somatic mutation
	flank_mutation_prob("lr_result.Rdata", "inner_anno_chr1.tsv", "outer_anno_chr1.tsv", "inner_pos_prob.tsv", "outer_pos_prob.tsv")

5. Run clustered_mut_prob() to combine the inner mutations and outer mutations for each recurrent mutation
	clustered_mut_prob("inner_pos_prob.tsv", "outer_pos_prob.tsv", "inner_mut_chr1.bed", "outer_mut_chr1.bed", "list_result.Rdata")
	
6. Run cal_clustered_score() to calculate clustered score in noncoding regions
	cal_clustered_score("list_result.Rdata")
	
## Getting Help ##
You mail the author(zhangjing1@shanghaitech.edu.cn) if you have any problem 
	
	
	
	