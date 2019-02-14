
#Example Context: Simple Linear Regression 
#To simplify things, we will sample for the residual variance and use theb

needed_packages = c("invgamma", "mvtnorm", "coda", "truncnorm")
for(i in 1:length(needed_packages)){
  haspackage = require(needed_packages[i], character.only = TRUE)
  if(haspackage == FALSE){
    install.packages(needed_packages[i])
  }
  library(needed_packages[i], character.only = TRUE)
}

# import data ==============================================================================================

setwd(dir = "Lectures/08-MCMC-Algorithms/")

# so I will use data from our previous example:
DietData = read.csv(file = "DietData.csv")

FullModel = lm(formula = WeightLB ~ HeightIN + factor(DietGroup) + HeightIN:factor(DietGroup), data=DietData)
summary(FullModel)

# Algorithm setup =======================================================================================================

nBeta = nrow(summary(FullModel)$coef)

nIterations = 10000
nChains = 4
nBurnin = 5000

# for tuning the chain (adapting):
maxTuneChain = 10
nTuneIterations = 500
tuneTarget = c(.40,.55)

# Gibbs Sampling for Betas ==============================================================================================

#set priors: normal for intercept/slope; inverse gamma for residual variance
betaMean.0 = 0
betaVariance.0 = 1000000 # a big variance for beta is needed as all will share a common prior mean

#put into matrices for use in calculations below
betaMeanVec.0 = matrix(data = rep(betaMean.0, nrow(summary(FullModel)$coef)), ncol = 1)
betaCov.0 = diag(nBeta)*betaVariance.0
betaInvCov.0 = solve(betaCov.0)

# Metropolis-Hastings Sampling for Sigma2 ===============================================================================

set.seed(seed = 12022019)

# We will draw Sigma2 from a proposal distribution that is truncated normal with mean = current sigma2 and variance that is adapted

# for sigma2_e (we will not use tau this time)
sigma2.0 = var(DietData$WeightLB)
nu.0 = 1

alpha.sigma2.0 = nu.0/2
beta.sigma2.0 = (nu.0*sigma2.0)/2

# rename data into something more useable
Y = matrix(data = DietData$WeightLB, ncol=1)
X = model.matrix(FullModel)
N = nrow(X)

# create chain list object
posterior = list()

# for progress bar
iterationCount = 0 

totalIterations =  nChains*(maxTuneChain*nTuneIterations + nIterations)
pb = txtProgressBar(min = 0, max = 1)

propVarSummary = NULL

chainNum = 1
for (chainNum in 1:nChains){
  cat("\n")
  
  # set a placeholder for the chain values
  chain = matrix(data = NA, nrow = nIterations, ncol = ncol(X)+1)
  colnames(chain) = c(colnames(X), "sigma2.e")
  
  #set starting values--randomly from values of the prior for now but can be changed
  Beta = matrix(data = rnorm(n = nBeta, mean = betaMean.0, sd = sqrt(betaVariance.0)), nrow = nBeta, ncol = 1)
  sigma2 = rtruncnorm(n = 1, a = 0, b = Inf, mean = summary(FullModel)$sigma^2, sd = 100)

  # initial proposal variance
  sigma2.propVar = 100
  
  # for counting iterations
  tune = 0 
  while (tune <= maxTuneChain){
    
    if (tune > 0){
      
      # convert number of acceptances to proportion (+1 in denominator to ensure qnorm() has no errors)
      propAccept = propAccept/(nTuneIterations+1)
      
      # see if we need to change the size of the variance of the proposal distribution
      if (propAccept < tuneTarget[1] | propAccept > tuneTarget[2] | tune == maxTuneChain){
        propVarSummary = rbind(propVarSummary, data.frame(chain = chainNum, tune = tune, propVar = sigma2.propVar, propAccept = propAccept))
        # need to adjust variance (using methods described by https://support.sas.com/documentation/cdl/en/statug/63962/HTML/default/viewer.htm#statug_mcmc_sect024.htm)
        sigma2.propVar = sigma2.propVar*qnorm(p = (mean(tuneTarget)/2))/qnorm(p = (propAccept/2))
        
        # reset counter using 1 so will never lead to an error in qnorm()
        propAccept = 1
        currentState = "tuning"
        
        phaseIterations = nTuneIterations
        
        tune = tune + 1
      } else {
        # variance is great...start sampling chain
        
        propVarSummary = rbind(propVarSummary, data.frame(chain = chainNum, tune = tune, propVar = sigma2.propVar, propAccept = propAccept))
        
        tune = maxTuneChain + 1
        currentState = "sampling"
        phaseIterations = nIterations
      }
      
    } else {
      # for first pass through algorithm
      
      tune = tune + 1
      propAccept = 0
      phaseIterations = nTuneIterations
    }
    
    iter = 1
    for (iter in 1:phaseIterations){
      
      #sample beta vector -- note: methods here are not best numerically and with respect to efficiency but are shown for didactic purposes
      BetaCov.n = solve(betaInvCov.0 + (1/sigma2)*t(X) %*% X)
      BetaMean.n = solve(betaInvCov.0 + (1/sigma2)*t(X) %*% X) %*% (betaInvCov.0 %*% betaMeanVec.0 + (1/sigma2)*t(X) %*% Y)
      
      #as return from rmvnorm is transposed from how we use beta in the above equations
      Beta = t(rmvnorm(n = 1, mean = BetaMean.n, sigma = BetaCov.n))
      
      #sample sigma2 using Metropolis-Hastings
      
      # propose new value of sigma2
      propSigma2 = rtruncnorm(n = 1, a = 0, b = Inf, mean = sigma2, sd = sqrt(sigma2.propVar))
      
      # calculate sample likelihood (using log)
      curLogLike = sum(dnorm(x = Y, mean = X %*% Beta, sd = sigma2, log = TRUE))
      propLogLike = sum(dnorm(x = Y, mean = X %*% Beta, sd = propSigma2, log = TRUE))
      
      # calculate prior log likelihoods
      curPriorLike = dinvgamma(x = sigma2, shape = alpha.sigma2.0, rate = beta.sigma2.0, log = TRUE)
      propPriorLike = dinvgamma(x = propSigma2, shape = alpha.sigma2.0, rate = beta.sigma2.0, log = TRUE)
      
      # calculate transition likelihoods
      qPropGivenCur = log(dtruncnorm(x = propSigma2, a = 0, b = Inf, mean = sigma2, sd = sqrt(sigma2.propVar)))
      qCurGivenProp = log(dtruncnorm(x = sigma2, a = 0, b = Inf, mean = propSigma2, sd = sqrt(sigma2.propVar)))
      
      # calculate log of MH ratio:
      logR = (propLogLike + propPriorLike + qCurGivenProp) - (curLogLike + curPriorLike + qPropGivenCur)
      
      # determine if we accept or reject
      if (log(runif(1)) < logR){
        # accept
        sigma2 = propSigma2
        propAccept = propAccept + 1
      } 
      
      # if during sampling phase, save parameters in chain
      if (tune == (maxTuneChain+1)) chain[iter,] = cbind(t(Beta), sigma2)
      currentIteration  =  ((tune-1)*(nTuneIterations) + iter) + (chainNum-1)*(maxTuneChain*nTuneIterations + nIterations)
      setTxtProgressBar(pb = pb, value = currentIteration/totalIterations)
    }
  }
  
  # add chain to chains list and convert to an mcmc object for easy coda processing
  posterior[[chainNum]] = as.mcmc(chain)
}
close(pb)

# finish by converting MCMC chain to MCMC.list 
posterior = mcmc.list(posterior)

#summarized chain:
plot(posterior)
summary(posterior)

# # trim off burnin and show chain
posterior = window(x = posterior, start = nBurnin+1, stop = nIterations)
plot(posterior)
summary(posterior)

propVarSummary
