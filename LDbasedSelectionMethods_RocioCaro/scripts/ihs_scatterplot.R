library("ggplot2")

YRI <- read.delim("YRI.ihs.out.100bins.norm", header=FALSE)

percentiles <- quantile(abs(YRI$V7), c(0.95, 0.99))

ggplot(YRI, aes(x=V2, y=V7)) + 
  geom_point(size=1) + 
  geom_hline(aes(yintercept=percentiles[1]), colour="blue") +
  annotate("text", x= Inf, y = percentiles[1], label="5%*", vjust=-0.5, hjust=1) +
  geom_hline(aes(yintercept=percentiles[2]), colour="red") +
  annotate("text", x= Inf, y = percentiles[2], label="1%*", vjust=-0.5, hjust=1) +
  geom_hline(aes(yintercept=-percentiles[1]), colour="blue") +
  annotate("text", x= Inf, y = -percentiles[1], label="5%*", vjust=-0.5, hjust=1) +
  geom_hline(aes(yintercept=-percentiles[2]), colour="red") +
  annotate("text", x= Inf, y = -percentiles[2], label="1%*", vjust=-0.5, hjust=1) +
  theme(legend.position="none") + xlab("Position") + ylab("iHS") + ggtitle("iHS for YRI")
  
ggsave("YRI_scatterplot.pdf", width = 10, height = 4)


CEU <- read.delim("CEU.ihs.out.100bins.norm", header=FALSE)

percentiles <- quantile(abs(CEU$V7), c(0.95, 0.99))

ggplot(CEU, aes(x=V2, y=V7)) + 
  geom_point(size=1) + 
  geom_hline(aes(yintercept=percentiles[1]), colour="blue") +
  annotate("text", x= Inf, y = percentiles[1], label="5%*", vjust=-0.5, hjust=1) +
  geom_hline(aes(yintercept=percentiles[2]), colour="red") +
  annotate("text", x= Inf, y = percentiles[2], label="1%*", vjust=-0.5, hjust=1) +
  geom_hline(aes(yintercept=-percentiles[1]), colour="blue") +
  annotate("text", x= Inf, y = -percentiles[1], label="5%*", vjust=-0.5, hjust=1) +
  geom_hline(aes(yintercept=-percentiles[2]), colour="red") +
  annotate("text", x= Inf, y = -percentiles[2], label="1%*", vjust=-0.5, hjust=1) +
  theme(legend.position="none") + xlab("Position") + ylab("iHS") + ggtitle("iHS for CEU")
  
ggsave("CEU_scatterplot.pdf", width = 10, height = 4)



