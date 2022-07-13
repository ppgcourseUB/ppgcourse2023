#This script helps you to fit the priors for theta and tau in BPP 4.0
#according to the distribution of your observed maximum theta per locus and pairwise divergence parameters.
#Take into account that the priors for theta and tau in BPP 4.0 follow a inverse gamma distribution (read BPP tutorial).
#First, you should have already estimated  "maximum theta per locus" and "pairwise divergence" by means the "3_EDCchopped_locifile.r" script

#Get the .csv file from the following "3_EDCchopped_locifile.r" script
stat.file <- read.csv('stat_chopped_loci.csv');

##THETA
par(mfrow = c(2,1));
hist(stat.file[,2], main = 'Maximum Theta per locus', xlab = 'Theta value (intraspecific)', xlim = c(-0.005,0.08), breaks = c(seq(-0.005, 0.8, by = 0.005)));
abline(v = quantile(stat.file[,2], 0.025), lwd = 1.5, lty = 2, col = 'red');
abline(v = mean(stat.file[,2], 0.5), lwd = 2, col = 'blue');
abline(v = quantile(stat.file[,2], 0.975), lwd = 1.5, lty = 2, col = 'red');

#Calculate the mean 2.5 and 97.5 quantiles of the distribution of the observed theta parameter
cat('The mean value is', mean(stat.file[,2]))
cat('The 2.5% and 97.5% quantile values are', quantile(stat.file[,2], 0.025), 'and', quantile(stat.file[,2], 0.975))


#Prior setting for theta...

library("invgamma")
a=3; b=0.04; #Change the value of b parameter. Keep alpha = 3 following recommendations from BPP tutorial
curve(dinvgamma(x,a,b), from=0, to=0.08)
abline(v = qinvgamma(c(0.5),a,b), lwd = 1.5, lty = 2, col = 'blue');
abline(v = qinvgamma(c(0.025),a,b), lwd = 1.5, lty = 2, col = 'red');
abline(v = qinvgamma(c(0.975),a,b), lwd = 1.5, lty = 2, col = 'red');

cat('The mean value is', qinvgamma(c(0.5),a,b)); #Calculate the mean value for a given invgamma distribution
cat('The 2.5% and 97.5% quantile values are', qinvgamma(c(0.025,0.975),a,b)) #Calculate the 2.5 and 97.5 quantiles for a given invgamma distribution

#Modify the b parameter until you fit the invgamma distribution to your observed data distribution
#Check your the the 2.5, 97.5 quantiles of your observed data distribution. Those values should be included within
#the 2.5 and 97.5 quantiles values for a the above invgamma distribution


##TAU
par(mfrow = c(2,1));
hist(stat.file[,3]*100, main = 'Maximum p-distance', xlab = '% pairwise divergence (interspecific)', xlim = c(-1,15), breaks = c(seq(-1, 20, by = 1)));
abline(v = quantile(stat.file[,3]*100, 0.025), lwd = 1.5, lty = 2, col = 'red');
abline(v = mean(stat.file[,3]*100, 0.5), lwd = 2, col = 'blue');
abline(v = quantile(stat.file[,3]*100, 0.975), lwd = 1.5, lty = 2, col = 'red');

par(mfrow = c(2,1));
hist(stat.file[,3], main = 'Maximum p-distance', xlab = 'Pairwise divergence (interspecific)', xlim = c(-0.01,0.15), breaks = c(seq(-0.01, 0.20, by = 0.01)));
abline(v = quantile(stat.file[,3], 0.025), lwd = 1.5, lty = 2, col = 'red');
abline(v = mean(stat.file[,3], 0.5), lwd = 2, col = 'blue');
abline(v = quantile(stat.file[,3], 0.975), lwd = 1.5, lty = 2, col = 'red');


#Calculate the mean value, abd 2.5 and 97.5 quantiles of the distribution of the observed theta parameter
cat('The mean value is (%)', mean(stat.file[,3]*100));
cat('The 2.5% and 97.5% quantile values are (%)', quantile(stat.file[,3]*100, 0.025), 'and', quantile(stat.file[,3]*100, 0.975))

cat('The mean value is', mean(stat.file[,3]));
cat('The 2.5% and 97.5% quantile values are', quantile(stat.file[,3], 0.025), 'and', quantile(stat.file[,3], 0.975))

#Prior setting for tau...

library("invgamma")
a=3; b=0.08; #Change the value of b parameter. Keep alpha = 3 following recommendations from BPP tutorial
curve(dinvgamma(x,a,b), from=0, to=0.15)
abline(v = qinvgamma(c(0.5),a,b), lwd = 1.5, lty = 2, col = 'blue');
abline(v = qinvgamma(c(0.025),a,b), lwd = 1.5, lty = 2, col = 'red');
abline(v = qinvgamma(c(0.975),a,b), lwd = 1.5, lty = 2, col = 'red');

cat('The mean value is', qinvgamma(c(0.5),a,b)) #Calculate the mnean value for a given invgamma distribution
cat('The 2.5% and 97.5% quantile values are', qinvgamma(c(0.025,0.975),a,b)) #Calculate the 2.5 and 97.5 quantiles for a given invgamma distribution

#Modify the b parameter until you fit the invgamma distribution to your observed data distribution
#Check your the the 2.5, 97.5 quantiles of your observed data distribution. Those values should be included within
#the 2.5 and 97.5 quantiles values for a the above invgamma distribution