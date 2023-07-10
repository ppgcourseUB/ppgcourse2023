#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir "00_Gene_filtering"
path_to_file <- getActiveDocumentContext()$path
script_wd <- paste( dirname( path_to_file ), "/", sep = "" )
wd <- gsub( pattern = "/scripts", replacement = "", x = script_wd )
setwd( wd )

#--------------#
# LOAD OBJECTS #
#--------------#
# Load tree 
raw_tt <- ape::read.tree( file = "../../00_data_formatting/00_raw_data/cytb_rooted_bl.tree" )

#-------#
# TASKS #
#-------#
# 1. Find tree height. You can use the function `phytools::nodeHeights` to
#    calculate all the heights of the tree. Then, we can extract the maximum
#    height calculated, which will correspond to the length from the root to 
#    the highest tip.
tree_height <- max( phytools::nodeHeights( raw_tt ) ) # 0.3801525

# 2. Get the mean of the calibration set for the root to have the 
#    time for the speciation event at the root, what we will use 
#    to estimate the mean evolutionary rate later. The first set of calibrations
#    has a soft-bound calibration to constraint the root age (i.e., a uniform
#    distribution with a minimum of 37.71 Ma and a 
#    maximum of 66.09 Ma with soft bounds). The average in time unit = 100 Ma 
#    is then 0.519 * 100 Ma. The mean root age following the second set of 
#    calibrations is 0.2338 * 100 Ma (min. age = 20.16 Ma | max. age = 26.6 Ma).
root_age_cal1 <- mean( c(0.3771,0.6609) ) # 0.519 * 100 Ma
root_age_cal2 <- mean( c(0.2016,0.266) )  # 0.2338 * 100 Ma 
# 3. Estimate mean rate based on the two different time units
#
#    tree_height = mean_rate x root_age --> mean_rate = tree_height / root_age
mean_rate_cal1 <- tree_height / root_age_cal1 # 0.7324711 subst/site per time unit
mean_rate_cal2 <- tree_height / root_age_cal2 # 1.625973 subst/site per time unit
# If we want to know the mean rate in subst/site/year, we apply the time unit. We
# We should get the same estimate regardless the time unit used:
#
# Calib1: Time unit 100 May (10^8y): 0.7324711 subst/site/10^8 = 
#                                    = 7.32e-9 subst/site/year
# Cali2: Time unit 100 May (10^8y): 1.625973 subst/site/10^8 = 
#                                    = 1.63e-9 subst/site/year

# 4. Now, we can build the gamma distribution given that we now have an estimate
#    for the mean rate. We will also use `alpha = 2` as we will start with a 
#    vague distribution. Nevertheless, if you were very sure about the mean 
#    rate, you could build a more constraint prior.
#
#    mean_Gamma = mean_rate = alpha / beta --> beta = alpha / mean_rate
alpha <- 2
beta_cal1  <- alpha/mean_rate_cal1 # 2.730483 ~ 2.7
beta_cal2  <- alpha/mean_rate_cal2 # 1.230033 ~ 1.2

# We can plot these distributions
curve( dgamma( x, shape = 2, rate = beta_cal1 ), from = 0, to = 3, col = "black" )
legend( "topright", legend = c( "G(2,2.7) " ), 
        col = "black", lty = 1, box.lty = 2 )
curve( dgamma( x, shape = 2, rate = beta_cal2 ), from = 0, to = 8, col = "black" )
legend( "topright", legend = c( "G(2,1.2) " ), 
        col = "black", lty = 1, box.lty = 2 )

# 5. Plot the gamma distributions
if ( ! dir.exists( "out_RData" ) ){
  dir.create( "out_RData" )
}
pdf( file = "out_RData/gamma_dists.pdf", paper = "a4" )
par( mfrow = c(1,2))
curve( dgamma( x, shape = 2, rate = beta_cal1 ), from = 0, to = 8, col = "black" )
legend( "topright", legend = c( "G(2,2.7) " ), 
        col = "black", lty = 1, box.lty = 2 )
curve( dgamma( x, shape = 2, rate = beta_cal2 ), from = 0, to = 8, col = "black" )
legend( "topright", legend = c( "G(2,1.2) " ), 
        col = "black", lty = 1, box.lty = 2 )
dev.off()

