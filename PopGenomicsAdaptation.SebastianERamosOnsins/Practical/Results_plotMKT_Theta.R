source("./asymptoticMK_local.R")

############################################
# FUNCTIONS
############################################

#calculation of alpha from standard equation
calc.alpha <- function(Ps,Ds,Pn,Dn) {
  alpha <- 1 - (Ds/Dn)*(Pn/Ps)
  return(alpha)
} 
#calculation of different theta estimators
CalcThetaUnfolded <- function(sfs,w) {
  th <- 0
  for(i in 1:length(sfs)) {
    th <- th + w[i] * i * sfs[i]
  }
  th <- th/(sum(w))
  return(th)
}
#weights to estimate watterson estimate (based on the number of variants under SNM)
weight.watt.unfolded <-function(nsam) {
  w <- array(0,dim=c(floor(nsam-1)))
  for(i in 1:length(w)) {
    w[i] <- 1/i
  }
  w
}
#weights to estimate tajima (nucleotide diversity, PI) estimate
weight.taj.unfolded <-function(nsam) {
  w <- array(0,dim=c(floor(nsam-1)))
  for(i in 1:length(w)) {
    w[i] <- nsam-i
  }
  w
}
#weights to estimate fu&li estimate (singletons)
weight.fuli.unfolded <-function(nsam) {
  w <- array(0,dim=c(floor(nsam-1)))
  w[1] <- 1
  w
}
#weights to estimate fay&wu estimate (based on high frequencies)
weight.fw.unfolded <- function(nsam) {
  w <- array(0,dim=c(floor(nsam-1)))
  for(i in 1:length(w)) {
    w[i] <- i
  }
  w
}

############################################
# DEFINITIONS
############################################

pdf("Plots_scenarios_alpha_Theta.pdf")
#define sample and size of sfs (to avoid zeros)
nsam <- 25*2
#read files from slim output
slim_files <- system("ls *_slim_SFS_*.txt",intern=T)
#Define data frame to keep results from variability estimates
Theta.results <- array(0,dim=c(length(slim_files),9))
colnames(Theta.results) <- c("scenario","Theta.FuLi.Nsyn","Theta.FuLi.Syn",
                             "Theta.Watt.Nsyn","Theta.Watt.Syn",
                             "Theta.Taji.Nsyn","Theta.Taji.Syn",
                             "Theta.FayWu.Nsyn","Theta.FayWu.Syn")
Theta.results <- as.data.frame(Theta.results)

Alpha.theta.results <- array(0,dim=c(length(slim_files),6))
colnames(Alpha.theta.results) <- c("scenario","alpha.real","alpha.mkt","alpha.assym","alphaAs.CIL","alphaAs.CIH")
Alpha.theta.results <- as.data.frame(Alpha.theta.results)

############################
#RUN
############################
i <- 1
#f <- slim_files[i]
for(f in slim_files) {
  dat.sfs <- read.table(file=f,header=T,row.names=1)
  #calc div
  divergence <- data.frame(mi=dat.sfs[1,nsam+2],Di=dat.sfs[1,nsam+1],m0=dat.sfs[2,nsam+2],D0=dat.sfs[2,nsam+1])

  #############################################
  #Estimation of alpha from theta estimators
  #############################################
  #Theta estimators are LIKE summary of sections of the SFS: 
  #It is useful for non-massive datasets
  #Fu&Li based on singletons
  #Watterson based on variants (mainly on variants at lower frequency)
  #Tajima weighting more on intermediate frequencies
  #Fay and Wu weighting more on higher frequencies
  
  #weights
  w.fuli <- weight.fuli.unfolded(nsam)
  w.watt <- weight.watt.unfolded(nsam)
  w.taj  <- weight.taj.unfolded(nsam)
  w.fayw <- weight.fw.unfolded(nsam)
  
  #Different Estimates of Variability;
  Theta.FuLi.Nsyn  <- CalcThetaUnfolded(sfs=as.numeric(dat.sfs[1,1:(nsam-1)]),w=w.fuli)
  Theta.FuLi.Syn   <- CalcThetaUnfolded(sfs=as.numeric(dat.sfs[2,1:(nsam-1)]),w=w.fuli)
  Theta.Watt.Nsyn  <- CalcThetaUnfolded(sfs=as.numeric(dat.sfs[1,1:(nsam-1)]),w=w.watt)
  Theta.Watt.Syn   <- CalcThetaUnfolded(sfs=as.numeric(dat.sfs[2,1:(nsam-1)]),w=w.watt)
  Theta.Taji.Nsyn  <- CalcThetaUnfolded(sfs=as.numeric(dat.sfs[1,1:(nsam-1)]),w=w.taj)
  Theta.Taji.Syn   <- CalcThetaUnfolded(sfs=as.numeric(dat.sfs[2,1:(nsam-1)]),w=w.taj)
  Theta.FayWu.Nsyn <- CalcThetaUnfolded(sfs=as.numeric(dat.sfs[1,1:(nsam-1)]),w=w.fayw)
  Theta.FayWu.Syn  <- CalcThetaUnfolded(sfs=as.numeric(dat.sfs[2,1:(nsam-1)]),w=w.fayw)
  
  #Estimates of standard alpha based on different estimates of variability
  alpha.mkt.FuLi <- calc.alpha(Ps=Theta.FuLi.Syn,Ds=divergence[1,4],Pn=Theta.FuLi.Nsyn,Dn=divergence[1,2])
  alpha.mkt.Watt <- calc.alpha(Ps=Theta.Watt.Syn,Ds=divergence[1,4],Pn=Theta.Watt.Nsyn,Dn=divergence[1,2])
  alpha.mkt.Taji <- calc.alpha(Ps=Theta.Taji.Syn,Ds=divergence[1,4],Pn=Theta.Taji.Nsyn,Dn=divergence[1,2])
  alpha.mkt.FayWu <- calc.alpha(Ps=Theta.FayWu.Syn,Ds=divergence[1,4],Pn=Theta.FayWu.Nsyn,Dn=divergence[1,2])
  
  #Assymptotic estimation APPROACH
  daf.theta <- array(0,dim=c(4,3))
  colnames(daf.theta)  <- c("daf","Pi","P0")
  #the variability estimates substitute the values on daf
  daf.theta[,2] <- sapply(c(Theta.FuLi.Nsyn,Theta.Watt.Nsyn,Theta.Taji.Nsyn,Theta.FayWu.Nsyn),round,4)
  daf.theta[,3] <- sapply(c(Theta.FuLi.Syn, Theta.Watt.Syn, Theta.Taji.Syn, Theta.FayWu.Syn ),round,4)
  daf.theta <- as.data.frame(daf.theta)
  #estimate mean frequency values for each estimate...
  daf.theta[1,1] <- 1/sum(w.fuli) * sum(w.fuli*c(1:(nsam-1))) * 1/nsam
  daf.theta[2,1] <- 1/sum(w.watt) * sum(w.watt*c(1:(nsam-1))) * 1/nsam
  daf.theta[3,1] <- 1/sum(w.taj)  * sum(w.taj *c(1:(nsam-1))) * 1/nsam
  daf.theta[4,1] <- 1/sum(w.fayw) * sum(w.fayw*c(1:(nsam-1))) * 1/nsam
  aa <- NULL
  tryCatch(
    {
      aa <- asymptoticMK(d0=divergence$D0, d=divergence$Di, xlow=0, xhigh=1, df=daf.theta, true_alpha=NA, output="table")
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
  alpha.mkt.daf.theta <- c(alpha.mkt.FuLi,alpha.mkt.Watt,alpha.mkt.Taji,alpha.mkt.FayWu)
  
  #PLOT results
  plot(x=daf.theta[,1],y=alpha.mkt.daf.theta,type="p",xlim=c(0,1),ylim=c(min(-1,alpha.mkt.daf.theta),1),
       main=sprintf("ALPHA from Theta: %s \nTrue=%.3f MKTa=%.3f MKT=%.3f",f,true.alpha,aa$alpha_asymptotic,aa$alpha_original),
       xlab="Freq",ylab="alpha")
  abline(h=0,col="grey")
  abline(v=c(0,1),col="grey")
  abline(h=true.alpha,col="blue")
  if(aa$model=="linear") abline(a=aa$a,b=aa$b,col="red")
  if(aa$model=="exponential") {
    x=seq(1:19)/20
    lines(x=x,y=aa$a+aa$b*exp(-aa$c*x),type="l",col="red")
  }
  
  Theta.results$scenario[i] <- f
  Theta.results$Theta.FuLi.Nsyn[i]  <- Theta.FuLi.Nsyn
  Theta.results$Theta.FuLi.Syn[i]   <- Theta.FuLi.Syn
  Theta.results$Theta.Watt.Nsyn[i]  <- Theta.Watt.Nsyn
  Theta.results$Theta.Watt.Syn[i]   <- Theta.Watt.Syn
  Theta.results$Theta.Taji.Nsyn[i]  <- Theta.Taji.Nsyn
  Theta.results$Theta.Taji.Syn[i]   <- Theta.Taji.Syn
  Theta.results$Theta.FayWu.Nsyn[i] <- Theta.FayWu.Nsyn
  Theta.results$Theta.FayWu.Syn[i]  <- Theta.FayWu.Syn
  
  Alpha.theta.results$scenario[i] <- f
  Alpha.theta.results$alpha.real[i] <- true.alpha
  Alpha.theta.results$alpha.mkt[i] <- aa$alpha_original
  Alpha.theta.results$alpha.assym[i] <- aa$alpha_asymptotic
  Alpha.theta.results$alphaAs.CIL[i] <- aa$CI_low
  Alpha.theta.results$alphaAs.CIH[i] <- aa$CI_high
  
  i <- i + 1
}
dev.off()
write.table(x=Theta.results,file="Theta.results.txt",row.names=F,quote=F)
write.table(x=Alpha.theta.results,file="Alpha.theta.results.txt",row.names=F,quote=F)
