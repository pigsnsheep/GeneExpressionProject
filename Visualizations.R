ggplot(data=kept_standard_deviations, aes(x=x)) + geom_histogram(bins = 200) + geom_vline(xintercept = 0.7)
ggplot(data=genedata_means, aes(x=x)) + geom_histogram(bins = 200) + geom_vline(xintercept = 0.5)
