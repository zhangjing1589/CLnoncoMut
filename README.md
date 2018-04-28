# CLnoncoMut
The R package to prioritize functional noncoding mutations

### install ###
> library(devtools)
> install_github("zhangjing1589/CLnoncoMut")


### Running CLnoncoMut ###

1. Use function get_lr_model() to train logistic regression model of annotated genomic regions
   get_lr_model("all_annotation_chr1.tsv","all_mut_chr1.bed","lr_result.Rdata") 
   the file "all_annotation_chr1.tsv" contains genomic base pair features: base pair(bp), replicating time(reptime), TFBS signal(tfbs), GERP RS score(gerp), DNaseI signal(dnase), CpG Island(cpg), transcript status(gencode_anno) 
   the file "all_mut_chr1.bed" contains somatic mutations(single base substitution) in ICGC database 
   the file "lr_result.Rdata" is trainning result used to predict mutation probability for each genomic position

2. Use function flank_mutation_prob() to predict mutation probility for each mutation
   flank_mutation_prob("","","","","")
