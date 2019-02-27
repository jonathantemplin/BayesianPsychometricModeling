
#Example Context: Simple Linear Regression 
#To simplify things, we will sample for the residual variance and use theb

needed_packages = c("invgamma", "mvtnorm", "coda")
for(i in 1:length(needed_packages)){
  haspackage = require(needed_packages[i], character.only = TRUE)
  if(haspackage == FALSE){
    install.packages(needed_packages[i])
  }
  library(needed_packages[i], character.only = TRUE)
}

#Simulate data ==============================================================================================

setwd(dir = "Lectures/08-MCMC-Algorithms/")

# so I will use data from our previous example:
DietData = read.csv(file = "DietData.csv")

FullModel = lm(formula = WeightLB ~ HeightIN + factor(DietGroup) + HeightIN:factor(DietGroup), data=DietData)
summary(FullModel)

# Gibbs Sampling ==============================================================================================

nBeta = nrow(summary(FullModel)$coef)

nIterations = 10000
nChains = 4
nBurnin = 5000

#set priors: normal for intercept/slope; inverse gamma for residual variance
betaMean.0 = 0
betaVariance.0 = 1000000 # a big variance for beta is needed as all will share a common prior mean

#put into matrices for use in calculations below
betaMeanVec.0 = matrix(data = rep(betaMean.0, nrow(summary(FullModel)$coef)), ncol = 1)
betaCov.0 = diag(nBeta)*betaVariance.0
betaInvCov.0 = solve(betaCov.0)

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
chainNum = 1
pb = txtProgressBar(min = 0, max = 1)
totalIterations = nChains*nIterations

for (chainNum in 1:nChains){
  cat("\n")
  # set a placeholder for the chain values
  chain = matrix(data = NA, nrow = nIterations, ncol = ncol(X)+1)
  colnames(chain) = c(colnames(X), "sigma2.e")
  
  #set starting values--randomly from values of the prior for now but can be changed
  Beta = matrix(data = rnorm(n = nBeta, mean = betaMean.0, sd = sqrt(betaVariance.0)), nrow = nBeta, ncol = 1)
  sigma2 = rinvgamma(n = 1, shape = alpha.sigma2.0, rate = beta.sigma2.0)

  iter = 2
  for (iter in 1:nIterations){
    
    
    #sample beta vector -- note: methods here are not best numerically and with respect to efficiency but are shown for didactic purposes
    BetaCov.n = solve(betaInvCov.0 + (1/sigma2)*t(X) %*% X)
    BetaMean.n = solve(betaInvCov.0 + (1/sigma2)*t(X) %*% X) %*% (betaInvCov.0 %*% betaMeanVec.0 + (1/sigma2)*t(X) %*% Y)
    
    #as return from rmvnorm is transposed from how we use beta in the above equations
    Beta = t(rmvnorm(n = 1, mean = BetaMean.n, sigma = BetaCov.n))
    
    #sample sigma2
    alpha.sigma2.n = (N + alpha.sigma2.0)/2
    beta.sigma2.n = beta.sigma2.0 + (1/2)*t(Y-X %*% Beta) %*% (Y-X %*% Beta)
    sigma2 = rinvgamma(n = 1, shape = alpha.sigma2.n, rate = beta.sigma2.n)
    
    #save parameters in chain
    chain[iter,] = cbind(t(Beta), sigma2)
    setTxtProgressBar(pb = pb, value = (iter+(nChains-1)*nIterations)/totalIterations)
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

# trim off burnin and show chain
posterior = window(x = posterior, start = nBurnin+1, stop = nIterations)
plot(posterior)
summary(posterior)

