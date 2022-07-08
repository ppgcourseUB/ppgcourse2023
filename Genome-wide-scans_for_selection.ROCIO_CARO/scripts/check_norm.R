library("ggplot2")

YRI <- read.delim("YRI.ihs.out.100bins.norm", header=FALSE)

pdf("YRI_unstandardized.pdf")
ggplot(YRI) + geom_point(aes(x=V3, y=V6)) + ggtitle("Unstandardized iHS in YRI") + xlab("Derived allele frequency") + ylab("iHS")
dev.off()

pdf("YRI_standardized.pdf")
ggplot(YRI) + geom_point(aes(x=V3, y=V7)) + ggtitle("Standardized iHS in YRI") + xlab("Derived allele frequency") + ylab("iHS")
dev.off()


CEU <- read.delim("CEU.ihs.out.100bins.norm", header=FALSE)

pdf("CEU_unstandardized.pdf")
ggplot(CEU) + geom_point(aes(x=V3, y=V6)) + ggtitle("Unstandardized iHS in CEU") + xlab("Derived allele frequency") + ylab("iHS")
dev.off()

pdf("CEU_standardized.pdf")
ggplot(CEU) + geom_point(aes(x=V3, y=V7)) + ggtitle("Standardized iHS in CEU") + xlab("Derived allele frequency") + ylab("iHS")
dev.off()

