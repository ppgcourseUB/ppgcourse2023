library("ggplot2")

xpehh_rs_file <- read.delim(YRIvsCEU.xpehh.out.norm, header=FALSE)

percentiles <- quantile(abs(xpehh_rs_file$V4), c(0.95, 0.99))

ggplot(xpehh_rs_file[xpehh_rs_file$V2=="15",], aes(x=V3, y=V4)) + 
geom_point(size=1) + 
geom_hline(aes(yintercept=percentiles[1]), colour="blue") +
geom_hline(aes(yintercept=percentiles[2]), colour="red") +
geom_hline(aes(yintercept=-percentiles[1]), colour="blue") +
geom_hline(aes(yintercept=-percentiles[2]), colour="red") +
theme(legend.position="none") + xlab("position") + ylab("XPEHH")

ggsave("YRIvsCEU_xpehh_chr15_scatterplot.pdf", width = 10, height = 4)






