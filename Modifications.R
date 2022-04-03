stdev_cutoff <- 0.7
mean_cutoff <- 0.5

standard_cutoff_check <- genedata_stdevs > stdev_cutoff

mean_cutoff_check <- genedata_means > mean_cutoff

keep_gene_check <- (standard_cutoff_check | mean_cutoff_check)





