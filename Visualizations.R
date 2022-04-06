ggplot(data=drugdata_sds, aes(x=x)) + geom_histogram(bins = 200) + geom_vline(xintercept = 0.5)
ggplot(data=drugdata_means, aes(x=x)) + geom_histogram(bins = 200) + geom_vline(xintercept = 0.77)

