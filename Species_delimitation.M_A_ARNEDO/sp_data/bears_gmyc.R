##Load the package
library(splits)

###read the ML tree
bears_c1_MLtree<- read.tree (file="./sp_data/bears_c1_root.treefile")

##Make the tree ultrametric
##There are multiple ways to linearise your tree, either by using your ML preferred tree, or by using infernce methods that already incoporate time estimaiton, eg. BEAST
##below is an example using the fuction chronos in R package APE (see http://phylobotanist.blogspot.com/2018/04/time-calibrated-or-at-least-ultrametric.html for some additional functions)
mycalibration <- makeChronosCalib(bears_c1_MLtree, node="root", age.min=4.5, age.max=5.5) ##assume a root age of interval 4.5 to 5.5 ma
bears_c1_MLchrono <- chronos(bears_c1_MLtree, lambda =10, model = "correlated", calibration = mycalibration, control = chronos.control() )

##Confirm that your tree fullfills the gmyc requirements
is.ultrametric(bears_c1_MLchrono)
is.binary.tree(bears_c1_MLchrono)
is.rooted(bears_c1_MLchrono)

##Otherwise, you can use the following tools
##Remove identical sequences
#remove.terminal.zeros(bears_c1_MLchrono)
## Make tree fully bifurcating
#multi2di(bears_c1_MLchrono)

##define & remove outgroup
outgroup <- match(c("NC009970_Melursus_ursinus"),bears_c1_MLchrono$tip.label)
bears_c1_MLchrono_nout <-drop.tip(bears_c1_MLchrono,outgroup)

## Run the single threshold GMYC model
bears_gmyc <- gmyc(bears_c1_MLchrono_nout, method="single", interval=c(0, 100))

## Send results to a log file
sink("myfile.log", append=TRUE, split=TRUE)

## Summarize results
summary.gmyc(bears_gmyc)

## Close log file
sink()

## Generate a list with the group assignments
bears_gmyc_groups<-spec.list(bears_gmyc) 

## Estimate support
pdf("bears_gmyc0.pdf")
yule_support <- gmyc.support(bears_gmyc)       
is.na(yule_support[yule_support == 0]) <- TRUE # only show values for affected nodes
plot(bears_gmyc, cex=.6, no.margin=TRUE)          # plot the tree
nodelabels(round(yule_support, 2), cex=.7)     # plot the support values on the tree
dev.off()
