# Approximate Bayesian Computation (ABC) in practice

**Instructor:** Oscar Lao

### 1. Parameter estimation given a model

Given an observed number of heads and tails from a coin, determine the posterior probability of the frequency of heads.Assume we throw 50 times a coin, recording each time if we get head (1) or tail (0):

```Rrm(list=ls());
coin <- rep(0,50);
coin[sample(1:50,rbinom(1,50,0.2))] <- 1
```

> `coin` is a simulated dataset that we have generated in silico. 

A simulated dataset refers to a collection of observations that are produced by a recognized model, utilizing specific parameter values associated with that model. A model is a statistical model when it follows a probabilistic distribution. **In our case, can we already say which is the frequency of heads from that simulated coin?**

+ **Task 1**. In order to run an ABC, we need a simulator to generate simulated data under the model we are considering. In our case, we want to generate _n_ coin throws given a probability of being head

    ```R
    generate.data <- function(n, p)
    {
      coin <- rep(0,n);
      coin[sample(1:n,rbinom(1,n,p))] <- 1;
      return(coin);
    }
    ```
> n, is the number of coin throws we want to simulate. p, is the probability of being head

***

+ **Task 2**. Now we need to define the Summary Statistics we will use. Summary statistics must be informative (= must correlate) with the parameter we want to estimate. For example, we will define as a first summary statistic the variance of the dataset _d_, which contains the dataset from throwing a coin _n_ times:

    ```R
    estimate.sum.stat.variance <- function(d)
    {
      # compute the variance on the vector of observations d and return the value
      return(var(d));
    }
    ```

   **Can you think in another statistic that (maybe) is better for estimating the     probability of heads in a coin?**

    ```R
    estimate.sum.stat.your.proposal <- function(d)
    {
      # your code here
      return(mean(d));
    }
    ```

   Imagine we define another sum. stat. that does not correlate AT ALL with the  frequency of the parameter to estimate. For example. A value sampled at random from a normal distribution with mean 0 and standard deviation 10:

    ```R
    estimate.sum.stat.dummy <- function(d)
    {
      # your code here
      return(rnorm(1,0,10));
    }
    ```
***

+ **Task 3**. Now we need to implement a routine for generating new datasets and retrieving the summary statistics. For each parameter from our model, we will have to define a prior probability that will define our prior beliefs on the shape of the distribution of that parameter. In our case, we only have one parameter to estimate, the frequency of heads.

   * _Task 3.1_. Between which values can range the frequency of heads? Is this range limited between a minimum and a maximum, or can the minimum and maximum range between -infinite and infinite?
   
   * _Task 3.2_. Now we need to find a distribution that can model our prior beliefs. For example, if our prior belief about the frequency of heads in our coin is completely unknown, we can use a probability distribution that has the same probability in all the range. The uniform distribution would be fine. Produce the value from an uniform distribution between 0 and 1
    
    ```R
    simulate.parameter.uniform <- function()
    {
    return(runif(1));
    }
    ```
      
   Imagine we have prior information about the frequency of heads in the coin. In other words: we know that there is a value for the frequency of heads that is a priori more likely to be the one that applies to our coin. Assume that we believe that the most likely value is 0.5.  **Can you propose a prior distribution that ranges between the minimum and maximum of the parameter to estimate with a mode at 0.5?**

   > HINT: Check the shape of the beta distribution

   Produce the value from an uniform distribution between 0 and 1:
    
    ```R
    simulate.parameter.beta <- function(YOUR.A, YOUR.B)
    {
      return(rbeta(1, YOUR.A, YOUR.B));
    }
    ```
***

+ **Task 4**. Now we have to generate simulated datasets sampling the value of the parameter from the prior distribution. We will store the value of the parameter we have used from the prior and the summary statistics we compute from the generated dataset

   * Observed data
   
     ```R
     coin
     ```
   * sum stats of coin
    
     ```R
     sums.stats.observed <-c(estimate.sum.stat.variance(coin),estimate.sum.stat.your.proposal(coin),estimate.sum.stat.dummy(coin));

     n.generated.datasets <- 1000000

     sim.param <- matrix(nrow=n.generated.datasets, ncol = 1, -1);

     # Here we store the sum stat variance, your proposal, and the dummy variable
     sim.sumstat <- matrix(nrow=n.generated.datasets, ncol = 3, -1);

     # Generate 1000000 simulated datasets. For each dataset, compute the sum stats
     for(sim in 1:n.generated.datasets)
     {
      frequency.of.head <- simulate.parameter.uniform();
      sim.param[sim,1] <- frequency.of.head; 
      sim.data <- generate.data(length(coin), frequency.of.head);
      sim.sumstat[sim,1] <- estimate.sum.stat.variance(sim.data);
      sim.sumstat[sim,2] <- estimate.sum.stat.your.proposal(sim.data); 
      sim.sumstat[sim,3] <- estimate.sum.stat.dummy(sim.data);
      }

     # Lets plot the first summary statistic vs the parameter...

     plot(sim.param[1:1000,1], sim.sumstat.scaled[1:1000,1])

     # Lets plot the second summary statistic vs the parameter...

     plot(sim.param[1:1000,1], sim.sumstat.scaled[1:1000,2])

     # Lets plot the third summary statistic vs the parameter...

     plot(sim.param[1:1000,1], sim.sumstat.scaled[1:1000,3])
     ```

   In other words, we have some summary statistics that do not correlate with the parameter we want to estimate, and another that perfectly correlate. In fact, the frequency estimated in a sample is a SUFFICIENT STATISTIC of the frequency of the model. This mean that we should remove non-informative statistics in order to get unbiased estimations of the posterior:

   We can intuitively expect that the parameter value of simulations whose summary statistics are similar to the ones computed in the observed data could also have generated our observed data.

   Compute the Euclidean distance between the sum stats of the observed dataset and the simulated ones

    ```R
    distance.observed.simulated <- rowSums((sums.stats.observed - sim.sumstat)^2)

   Identify the 1000 simulations closest to the summary statistics of the observed dataset

    ```R
    threshold <- quantile(distance.observed.simulated, 1000/n.generated.datasets)
    ```

   Now pick the values used to generate the closest simulations to the observed data.

    ```R
    # THIS IS SAMPLED FROM OUR POSTERIOR DISTRIBUTION!!!!!
    posterior.distribution <- sim.param[distance.observed.simulated <= threshold,1];
    ```

   Once we have a sample of the posterior distribution, we need to compute centrality and dispersion statistics. Also, we can look the shape of the distribution

    ```R
    plot(density(posterior.distribution), main = "My posterior density")
    ```

   The first thing we observe is that the distribution is around 0.2, which is fine (this is in fact the value we used for the observed simulated dataset). However, it is a bit weird the plato we observe. It does not look bell shaped. Despite showing a bell shape is problem dependent (in other words, it is perfectly fine not observing a bell shape), in this case we MUST expect a "bell shaped" distribution. Basically, because theory says that the posterior distribution of a binomial distribution can be  estimated using a conjugate prior beta distribution with shape uniform and the posterior distribution has the format _beta(y + 1, n - y + 1)_ where y is the number of heads and n is the total number of observations. In our case:

    ```R
    lines(density(rbeta(10000, sum(coin) +1, length(coin)-sum(coin)+1)), col = "red", lwd = 2)
    ```

   WHAT IS GOING ON??? Could it be that the shape of the posterior distribution depends on how informative are the summary statistics? Lets repeat the analysis but excluding the last summary statistic, that in fact was not informative at all

    ```R
    distance.observed.simulated <- rowSums((sums.stats.observed[-3] - sim.sumstat[,-3])^2)

    threshold <- quantile(distance.observed.simulated, 1000/n.generated.datasets)

    posterior.distribution <- sim.param[distance.observed.simulated <= threshold,1];

    plot(density(posterior.distribution), main = "My posterior density")

    lines(density(rbeta(10000, sum(coin) +1, length(coin)-sum(coin)+1)), col = "red", lwd = 2)
    ```

   WHOA! Now it looks like the expected posterior distribution using a classical Bayesian Framework! However... When we compute the distance between sum statistics, we are not taking into account that each summary statistic has a different magnitude:

    ```R
    vars <- apply(sim.sumstat,2,var)
    ```

   We see that the first statistic is ~ 14 times smaller than the second statistic, and 16006 smaller than the third one...What does that mean when computing a distance? How can we fix this situation?

    ```R
    merged.sums.stats <- rbind(sums.stats.observed,sim.sumstat)
    ```

   Scale everything. Each column will have variance 1. Now we can COMPARE statistics as they are in the same magnitude

    ```R
    merged.sums.stats <- scale(merged.sums.stats);

    sums.stats.observed <- merged.sums.stats[1,];

    sim.sumstat <- merged.sums.stats[-1,];

    distance.observed.simulated <- rowSums((sums.stats.observed[-3] - sim.sumstat[,-3])^2)

    threshold <- quantile(distance.observed.simulated, 1000/n.generated.datasets)

    posterior.distribution <- sim.param[distance.observed.simulated <= threshold,1];

    plot(density(posterior.distribution), main = "My posterior density")
    ```

   **So, key points to take into account: BEWARE YOUR SUMMARY STATISTICS!!!!!**

***

So far we have implemented the classical "rejection" algorithm of ABC. Which are the steps of this algorithm? Initially, the only way to do ABC was to generate your own code, as we have done in the previous exercise. Now lets going to redo the same thing but with the ABC package, which allows other types of ABC algorithms and conducts all the preprocessing.

```R
library("abc");

# check ?abc

target <- sums.stats.observed;
param <- sim.param;
sumstat <- sim.sumstat;
tol <- 1000/nrow(sumstat);

# The previous lines of code condensed in a single line...
do.abc <- abc(target,param,sumstat,tol, method = "rejection")

plot(density(do.abc$unadj.values))
lines(density(rbeta(10000, sum(coin) +1, length(coin)-sum(coin)+1)), col = "red", lwd = 2)
```

The rejection algorithm is a bit unfair, if you think a bit. We can argue that even if we are picking the closest 1000 simulations, there are some of these simulations that are closer than others to the observed data. We should give more weight to these simulations than to the ones that are more far away. This is done with a local linear regression.

```R
do.abc <- abc(target,param,sumstat,tol, method = "loclinear")

# Take into account that we are now taking the ADJUSTED VALUES, NOT THE RAW VALUES (THAT WOULD BE REJECTION ALGORITHM)
plot(density(do.abc$adj.values));
lines(density(rbeta(10000, sum(coin) +1, length(coin)-sum(coin)+1)), col = "red", lwd = 2)
```

In this simple case, there is no real need to do that, as we get more or less the same output), now from the sampled from the posterior distribution we can report some statistics:

```R
# CENTRALITY STATISTICS:

mean(do.abc$adj.values)
median(do.abc$adj.values)

# DISPERSION STATISTICS:

quantile(do.abc$adj.values, probs = c(0.025,0.975))

library(HDInterval)

# The HDI is the interval which contains the required mass such that all points within the interval
# have a higher probability density than points outside the interval.

Highest.Posterior.Density.Intervals <- hdi(do.abc$adj.values, 0.8)

summary(Highest.Posterior.Density.Intervals)
```

**To which extend are we doing it well? Is the prediction of the mean close enough to the real value given the summary statistics we are using given the error rate we are using?**

```R
threshold.tol <- c(.005,.01, 0.05) # What if I pick the closest 0.005, 0.01 or 0.05 simulated datasets?

# Since it takes a while, we will take only 10 samples
cv.res.reg <- cv4abc(data.frame(Na=param[,1]), sumstat, nval=10, tols=threshold.tol, method="loclinear");
```

**To which extend the centrality statistics are close to the real value?** Predicted mean vs value of the parameter that generated the data

```R
n.reps <- 10;
diff.mean.real <- rep(n.reps,0);

for(rep in 1:n.reps)
{
  observed.data <- sample(1:nrow(sumstat),1);
  target <- sumstat[observed.data,];
  param <- sim.param[-observed.data,];
  tol <- 1000/length(param);
  do.abc <- abc(target,param,sumstat[-observed.data,],tol, method = "loclinear")
  diff.mean.real[rep] <- sim.param[observed.data] - mean(do.abc$adj.values);
}
```

**How big is going to be the difference between the mean estimated and the real value?**

```R
mean(abs(diff.mean.real))
```


### 2. Model comparison using abc

```R
rm(list=ls())
library("moments")
```

In the Bayesian paradigm, models are also parameters. This means that we can use the Bayesian framework to estimate the posterior distribution of each model

+ **Task 1**. Lets assume we have this data:

    ```R
    obs.data <- c(7.1858181,11.0938148,6.4044882,4.5047751,11.2077658,12.5797382,2.6840912,1.0640095,8.4082828,19.0201338,6.6116233,6.6285363,9.9166719,16.9371605,7.5659874,12.4425144,15.3466282,8.6929810,6.9484457,12.2287024,16.8991237,5.4697094,14.8373391,13.3426023,9.9992937,11.8750643,11.9554763,16.5623567,14.1896326,14.0065793,17.7305306,2.4636234,11.6579442,7.0606132,7.9101501,10.4307592,8.1211206,9.3771855,7.9051378,6.6950816,0.1616765,13.1589058,15.2026316,7.2795010,10.1026212,6.2491684,8.5669489,17.0469774,1.7225283,4.7860519,4.4052375,15.9330358,15.4676232,8.9812559,-0.7203585,4.2658463,14.6183584,10.8202862,9.4582488,0.7096979,18.1179716,12.1857587,3.0861436,7.3134261,5.1837287,4.8280980,6.9045024,3.9521058,3.5558720,12.5663511,4.3865937,21.0853330,15.0197656,6.4542300,8.5937835,18.8173465,3.7786510,13.3952906,6.2211788,0.3992227,11.4882243,6.4386650,9.9606916,3.8053607,12.2169972,3.2897962,17.6093680,11.7251748,8.7080017,10.5723543,9.2146612,12.0574657,6.5920890,6.1412366,10.7575179,2.5630478,13.9615675,7.6326555,12.0438835,9.6972823);

    plot(density(obs.data), main = "My observed data")
    ```
   We have two different statistical models that can potentially explain this data:
   
   (i) Normal distribution with parameters mu and sd
   _mu_ can range between 2 and 20. We will use a uniform prior
   _sd_ can range between 1 and 10. We will use a uniform prior

   (ii) Gamma distribution with parameters shape and scale
   _shape_ can range between 1 and 20. We will use a uniform prior
   _scale_ between 0.1 and 2. We will use a uniform prior

   **Which is the posterior probability that model 1 (normal) has produced the dataset compared to model 2 (gamma)?**
   As before, we need a simulator for normal distribution and for gamma distribution:

    ```R
    # Normal distribution
    simulator.model.1 <- function(n, mu, sd)
    {
      return(rnorm(n,mean = mu, sd = sd));
    }

    simulator.model.2 <- function(n, shape, scale)
    { 
      return(rgamma(n, shape = shape, scale = scale));
    }
    ```

   We will need to generate summary statistics. We will put all of them in a single row.

    ```R
    # Compute some summary statistics

    sum.stats.to.compute <- function(v)
    {
     ss <- rep(0,8);
  
     # mean
  
     ss[1] = mean(v);
  
     # median
  
     ss[2] = median(v);
  
     # mode-like
  
     kk <- hist(v,plot=F);
  
     if(length(kk$mids[kk$density==max(kk$density)])==1)
     {
       ss[3] = kk$mids[kk$density==max(kk$density)];    
     }
     else
     {
       ss[3] = kk$mids[kk$density==max(kk$density)][1];    
     }
  

     # sd
  
     ss[4] = sd(v);
  
     # min
  
     ss[5] = min(v);
  
     # max
  
     ss[6] = max(v);
  
     # skewness
  
     ss[7] = skewness(v);
  
     # kurtosis
  
     ss[8] = kurtosis(v);
  
     return(ss);
    }
    ```
   Now we generate sims of each model according to their prior probability lets assume we believe our data has been generated with model 2 with 70% probability.

    ```R
    obs.sums.stats <- sum.stats.to.compute(obs.data);
    n <- length(obs.sums.stats);
    n_sims <- 1000000;
    models <- rep(NA, n_sims);
    # store the simulated statistics
    sim.sumstat <- matrix(nrow=n_sims, ncol = 8, NA);

    for(s in 1:n_sims)
    {
      # in probability of 0.7, pick model 2
      if(runif(1) <= 0.7)
      {
       shape <- runif(1,1,20);
       scale <- runif(1,0.1,2);
       sim.dataset <- simulator.model.2(n, shape, scale);
       models[s] <- 2; 
       sim.sumstat[s,] <- sum.stats.to.compute(sim.dataset);
      }
      else
      {
      # in probability of 0.3, pick model 1
      mean.sim <- runif(1,2,20);
      sd.sim <- runif(1,1,10);
      sim.dataset <- simulator.model.1(n, mean.sim, sd.sim);
      models[s] <- 1; 
      sim.sumstat[s,] <- sum.stats.to.compute(sim.dataset);    
      }
    }

    threshold <- 1000/n_sims;
    target <- obs.sums.stats;
    models <- as.numeric(models);
    ```

   Using the classical rejection method:
   
    ```R
    model.posterior.probability.model <- postpr(target = target, index = models, sumstat = sim.sumstat, tol = threshold, method = "rejection", corr = FALSE)

    summary(model.posterior.probability.model)
    ```

   Using the logistic approximation:
   
   ```R
    model.posterior.probability.model <- postpr(target = target, index = models, sumstat = sim.sumstat, tol = threshold, method = "mnlogistic", corr = FALSE)
    summary(model.posterior.probability.model)
    ```

   Looks like our data comes from a normal distribution (in fact, they were generated from a N(10,5)). Try to repeat the analysis but now setting corr = TRUE. **What is doing this parameter?**

***

   * _Task 1.2_. We have run the abc without paying so much attention to the summary statistics. Lets check how the different sum stats look like

       ```R
       pairs(sim.sumstat[1:1000,])

       cor(sim.sumstat)
       ```

      Looks like we have a substantial number of variables that are highly correlated. For example, var 1 and var 2 have a correlation of almost ONE. **HOW CAN IT BE?**. HAVING REDUNDANT summary statistics can be a potential problem. There are different ways to fix this problem, all of them based on machine learning techniques. From unsupervised machine learning, we could apply a PCA and use the first k generated PCs as summary statistics. These PCs will be by definition uncorrelated.

       ```R
       pca.approach <- rbind(obs.sums.stats, sim.sumstat);

       pca.res <- prcomp(pca.approach, scale = TRUE)

       plot(pca.res)
       ```

      Just by eye, we can say that with the first three PCs we have covered a substantial amount of the total variance our observed summary statistics correspond to the first line of the transformed PCs Do the pca with only the first three PCs.

       ```R  
       obs.sums.stats <- pca.res$x[1,1:3];

       # and the summary statistics:

       sim.sumstat <- pca.res$x[-1,1:3];

       threshold <- 1000/n_sims;
       target <- obs.sums.stats;
       models <- as.numeric(models);
       ```

      Using the classical rejection method:
   
       ```R
       model.posterior.probability.model <- postpr(target = target, index = models, sumstat = sim.sumstat, tol = threshold, method = "rejection", corr = FALSE)

       summary(model.posterior.probability.model)
       ```

      Using the logistic approximation:
   
      ```R
       model.posterior.probability.model <- postpr(target = target, index = models, sumstat = sim.sumstat, tol = threshold, method = "mnlogistic", corr = FALSE)
       summary(model.posterior.probability.model)
      ```

      Results are now quite different...**Key point:** Including non-informative statistics and highly correlated summary statistics can BIAS your results...There are different ways to remove redundant and non informative data. One possibility is to use an artificial neural network. We pass the summary statistics and try to make a prediction of the parameter/model. Then, we use this prediction as SUMMARY STATISTIC.

***

+ **Task 2**. **What happens if the model that generated the data is not included in the model comparison. For example, What happens if our priors are OUT of the values that generated the data?**

    ```R
    obs.sums.stats <- sum.stats.to.compute(obs.data);
    n <- length(obs.sums.stats);
    n_sims <- 1000000;
    models <- rep(NA,n_sims);
    # store the simulated statistics
    sim.sumstat <- matrix(nrow=n_sims, ncol = 8, NA);

    for(s in 1:n_sims)
    {
      # in probability of 0.7, pick model 2
      if(runif(1) <= 0.7)
      {
        shape <- runif(1,1,20);
        scale <- runif(1,0.1,2);
        sim.dataset <- simulator.model.2(n, shape, scale);
        # gamma will be model 2
        models[s] <- 2; 
        sim.sumstat[s,] <- sum.stats.to.compute(sim.dataset);
      }
      else
      {
    ```
    # In probability of 0.3, pick model 1. HOWEVER, THE PRIOR IS FAAAAR AWAY FROM THE REGION WHERE THE PARAMETERS FROM THE OBSERVED DATA WERE GENERATED. The prior distribution of the mean is REALLY far away from the REAL value we used for generating our data...
    
    mean.sim <- runif(1,100,200);
    sd.sim <- runif(1,1,10);
    sim.dataset <- simulator.model.1(n, mean.sim, sd.sim);
    # normal will be model 1
    models[s] <- 1; 
    sim.sumstat[s,] <- sum.stats.to.compute(sim.dataset);    
      }
    }
    ```
   
   Ignoring the fact that some of the summary statistics are highly correlated:

    ```R
    threshold <- 1000/n_sims;
    target <- obs.sums.stats;
    models <- as.numeric(models);
    ``` 
   Using the classical rejection method
   
    ```R
    model.posterior.probability.model <- postpr(target = target, index = models, sumstat = sim.sumstat, tol = threshold, method = "rejection", corr = FALSE)

    summary(model.posterior.probability.model)
    ```

   Now it turns out that the "best model" is 1...**Key message:** THE PRIOR DISTRIBUTIONS OF THE PARAMETERS AFFECTS THE POSTERIOR OF THE MODELS. **Key message**: Because or prior belief is model 1, we can be more tempted to accept this result as it is.

***

+ **Task 3**. **How much can we trust the posterior probabilities of the parameters of a model?** After deciding that the best model is 1 (the gamma distribution) even if the right one is the 2, we can estimate the posterior distributions of the two parameters (shape and scale) from this model:

    ```R
    obs.sums.stats <- sum.stats.to.compute(obs.data);
    n <- length(obs.sums.stats);
    n_sims <- 1000000;
    # we will store here the shape and scale
    sim.param <- matrix(nrow=n_sims, ncol = 2, -1);

    sim.sumstat <- matrix(nrow=n_sims, ncol = 8, -1);

    # Generate 1000000 simulated datasets. For each dataset, compute the sum stats
    for(sim in 1:n_sims)
    {
      shape <- runif(1,1,20);
      scale <- runif(1,0.1,2);
      sim.param[sim,1] <- shape;
      sim.param[sim,2] <- scale;
      sim.dataset <- simulator.model.2(n, shape, scale); 
      sim.sumstat[sim,] <- sum.stats.to.compute(sim.dataset);
    }

    # Estimate the posterior distributions of shape and scale

    tol <- 1000/nrow(sim.param);

    do.abc <- abc(obs.sums.stats,sim.param,sim.sumstat,tol, method = "loclinear")

    # shape

    plot(density(do.abc$adj.values[,1]))

    # scale

    plot(density(do.abc$adj.values[,2]))
    ```

  **Do you notice something weird in the estimated posterior distributions when using local linear regression?** Right, in principle, the posterior distribution should be within the margins of the prior. However, we can get NEGATIVE values for the shape and scale!!!!!. We can FORCE THEM TO BE IN THE PRIOR RANGE by rescaling the parameters before we run the abc:

    ```R
    # Min and max boundary of each parameter
    logit.bounds <- matrix(nrow = 2, ncol = 2, 0)
    logit.bounds[1,] <- range(sim.param[,1]);
    logit.bounds[2,] <- range(sim.param[,2]);


    do.abc <- abc(obs.sums.stats,sim.param,sim.sumstat,tol, method = "loclinear",      transf = "logit", logit.bounds = logit.bounds);
    ```

   Now you can see that the values are within the bounds of the prior distributions. If we use the classical rejection method, then there is no need to do the re-scaling: by definition we will get the posterior in the range of the prior.

   Once we have the posteriors, we could attempt to generate simulated datasets from this posterior. If the model is the one that generated our data, we can expect that the output will be similar to the observed dataset...

    ```R
    number.of.simulated.datasets <- 1000

    my.simulated.datasets.from.posterior <- matrix(nrow = number.of.simulated.datasets, ncol = length(obs.data), 0);

    n <- length(obs.data)

    for(s in 1:nrow(my.simulated.datasets.from.posterior))
    {
      # sample a row at random
      row.i <- sample(1:nrow(do.abc$adj.values),1);
      shape <- do.abc$adj.values[row.i,1];
      scale <- do.abc$adj.values[row.i,2];
      my.simulated.datasets.from.posterior[s,] <- simulator.model.2(n, shape, scale);
    }
    ```

   Now we have simulated data that should look like the observed data. A very quick and simple check out is to compare the simulated and the observed dataset. "Comparison" can be done by different ways. One would be using the SAME OR OTHER summary statistics. To make it simple, we will use the same summary statistics:

    ```R
    obs.sums.stats <- sum.stats.to.compute(obs.data);
    sim.sums.stats <- matrix(nrow = nrow(my.simulated.datasets.from.posterior), ncol = length(obs.sums.stats), 0);

    for(r in 1:nrow(sim.sums.stats))
    {
      sim.sums.stats[r,] <- sum.stats.to.compute(my.simulated.datasets.from.posterior[r,]);
    }

    merge.sim.sums <- cbind(obs.sums.stats,sim.sums.stats);

    pca.replicas <- prcomp(merge.sim.sums, scale = T);

    plot(pca.replicas$x[,1], pca.replicas$x[,2])

    labs <- rep("",nrow(merge.sim.sums));
    labs[1] <- "OBS";

    points(pca.replicas$x[1,1], pca.replicas$x[1,2], pch = 19, col = "red") text(pca.replicas$x[,1], pca.replicas$x[,2], labs)
    ```

   We can see that the observed data falls within the cloud, but not close to the center (0,0)...**KEY MESSAGE:** Check always the performance of the posterior on simulated data. At least with the summary statistics you use!!!! BEWARE YOUR PRIOR DISTRIBUTIONS!!!!!!

# 3. Special case with demographic models

Following the ABC paradigm, you will need:

1) A demographic simulator to generate genomic data.
    a. Frequency based data
    b. Structural variants
    c. Sequence based data
    d. ...
    Popular simulators: msprime, fsc,...

 2) Define demographic models, and parameters. Define prior distributions

 3) Define a set of summary statistics.
    a. Related to SFS
    b. Related to LD
    c. Versions based on deep learning can take as input the sequence
    d. ...

 4) Define the number of simulations. Depends on the computational complexity of models and computational resources. Can take DAYS.

If you want to implement these approaches for your problem, please contact us: oscar.lao@ibe.upf-csic.es









 
