source("./asymptoticMK_local.R")

############################################
# FUNCTIONS
############################################

#calculation of alpha from standard equation
calc.alpha <- function(Ps,Ds,Pn,Dn) {
  alpha <- 1 - (Ds/Dn)*(Pn/Ps)
  return(alpha)
} 
# D0,P0 are neutral, Di, Pi are functional
#define distribution of allele frequencies from smaller intervals (to avoid zeros)
calc2.daf <- function(tab.dat,nsam,intervals) {
  daf.red <- data.frame(daf=as.numeric(sprintf("%.3f",c(c(0:(intervals-1))/intervals,1))),Pi=0,P0=0)
  tab.dat <- rbind(tab.dat,c(0,0))
  for(i in c(1:intervals)) {
    daf.red$Pi[i] <- sum(tab.dat[(as.integer(daf.red$daf[i]*nsam)+1):(as.integer(daf.red$daf[i+1]*nsam)),1])
    daf.red$P0[i] <- sum(tab.dat[(as.integer(daf.red$daf[i]*nsam)+1):(as.integer(daf.red$daf[i+1]*nsam)),2])
  }
  daf.red <- daf.red[-(intervals+1),]
  return(daf.red)
}

############################################
# DEFINITIONS
############################################

pdf("Plots_scenarios_alpha.pdf")
#define sample and size of sfs (to avoid zeros)
nsam <- 25*2
intervals <- 50
#read files from slim output
slim_files <- system("ls *_slim_SFS_*.txt",intern=T)
#Define data frame to keep results from SFS
Alpha.results <- array(0,dim=c(length(slim_files),6))
colnames(Alpha.results) <- c("scenario","alpha.real","alpha.mkt","alpha.assym","alphaAs.CIL","alphaAs.CIH")
Alpha.results <- as.data.frame(Alpha.results)

############################
#RUN
############################
i <- 1
#f <- slim_files[i]
for(f in slim_files) {
  dat.sfs <- read.table(file=f,header=T,row.names=1)

  ############################
  # USING SFS
  ############################
  #calc daf and div
  tab.dat <- t(dat.sfs[,1:(nsam-1)])
  daf <- calc2.daf(tab.dat,nsam,intervals)
  #to avoid ZERO values, instead (or in addition) of doing less intervals, add 1 to freqs=0
  daf[daf[,2]==0,2] <- 1
  daf[daf[,3]==0,3] <- 1
  #calc div
  divergence <- data.frame(mi=dat.sfs[1,nsam+2],Di=dat.sfs[1,nsam+1],m0=dat.sfs[2,nsam+2],D0=dat.sfs[2,nsam+1])
  
  #estimate MKTa from SFS
  aa <- NULL
  tryCatch(
    {
      aa <- asymptoticMK(d0=divergence$D0, d=divergence$Di, xlow=0.1, xhigh=0.9, df=daf, true_alpha=NA, output="table")
    },
    error = function(e) {
      message(sprintf("Error calculating MKTa for observed data in file %s",f))
    }
  )
  aa
  aa$alpha_asymptotic
  aa$alpha_original
  true.alpha <- dat.sfs[1,nsam+3] / dat.sfs[1,nsam+1]
  true.alpha
  #calculate alpha for each frequency separately (assumed same nsam at Syn and Nsyn)
  alpha.mkt.daf <- calc.alpha(Ps=daf[,3],Ds=divergence[1,4],Pn=daf[,2],Dn=divergence[1,2])
  alpha.mkt.daf
  
  #Plot results
  plot(x=daf[,1],y=alpha.mkt.daf,type="p",pch=20,xlim=c(0,1),ylim=c(min(-1,alpha.mkt.daf),1),
       main=sprintf("ALPHA: %s \nTrue=%.3f MKTa=%.3f MKT=%.3f",f,true.alpha,aa$alpha_asymptotic,aa$alpha_original),
       xlab="Freq",ylab="alpha")
  abline(h=0,col="grey")
  abline(v=c(0.1,0.9),col="grey")
  abline(h=true.alpha,col="blue")
  if(aa$model=="linear") abline(a=aa$a,b=aa$b,col="red")
  if(aa$model=="exponential") {
    x=seq(1:19)/20
    lines(x=x,y=aa$a+aa$b*exp(-aa$c*x),type="l",col="red")
  }
  Alpha.results$scenario[i] <- f
  Alpha.results$alpha.real[i] <- true.alpha
  Alpha.results$alpha.mkt[i] <- aa$alpha_original
  Alpha.results$alpha.assym[i] <- aa$alpha_asymptotic
  Alpha.results$alphaAs.CIL[i] <- aa$CI_low
  Alpha.results$alphaAs.CIH[i] <- aa$CI_high
  
  i <- i + 1
}
dev.off()
write.table(x=Alpha.results,file="Alpha.results.txt",row.names=F,quote=F)
