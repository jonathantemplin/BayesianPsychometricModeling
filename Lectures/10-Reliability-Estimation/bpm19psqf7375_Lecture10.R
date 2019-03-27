needed_packages = c("rjags", "HDInterval", "mvtnorm")
for(i in 1:length(needed_packages)){
  haspackage = require(needed_packages[i], character.only = TRUE)
  if(haspackage == FALSE){
    install.packages(needed_packages[i])
  }
  library(needed_packages[i], character.only = TRUE)
}


set.seed = 1071
# 
nItems = 3
nObs = 1000
# 
# simulating data from a one-factor CFA model 
Mu = runif(n = nItems, min = 1, max = 3)
Lambda = matrix(data = runif(n = nItems, min = .5, max = 1.5), nrow = nItems, ncol = 1)
Phi = matrix(data = 1, nrow = 1, ncol = 1)
Psi = diag(runif(n = nItems, min = 0, max =1))

data = rmvnorm(n = nObs, mean = Mu, sigma = Lambda %*% Phi %*% t(Lambda) + Psi)

nchains = 4
niter = 5000
nburnin = 2000
nadapt = 2000
nthin = 1

# specification of prior values for measurement model parameters:
#   item means
mu.mean.0 = 3
mu.variance.0 = 1000000
mu.precision.0 = 1 / mu.variance.0

#   Factor loadings
lambda.mean.0 = 0
lambda.variance.0 = 1000000
lambda.precision.0 = 1 / lambda.variance.0

# unique variances
psi.df.0 = 1
psi.var.0 = 1
psi.alpha.0 = psi.df.0 / 2
psi.beta.0 = (psi.df.0 * psi.var.0) / 2

# values for prior for factor variance (based on variance of marker item)
factor.df.0 = 1
factor.var.0 = 1
factor.alpha.0 = factor.df.0/2
factor.beta.0 = (factor.df.0*factor.var.0)/2

# marker item:
model01.syntax = "
model{
  # measurement model specification
  for (person in 1:N){
    for (item in 1:I){
      mean[person, item] = mu[item] + lambda[item]*xfactor[person]
      X[person, item] ~ dnorm(mean[person,item], inv.psi[item])    
    }
  }

  # prior distributions for the factor:
  for (person in 1:N){
    xfactor[person] ~ dnorm(0, factor.precision)
  }

  # prior distributions for the measurement model mean/precision parameters
  for (item in 1:I){
    mu[item] ~ dnorm(mu.mean.0, mu.precision.0)
    inv.psi[item] ~ dgamma(psi.alpha.0, psi.beta.0)
  }

  # prior distributions for the loadings (except the first loading, which is fixed to 1.0)
  lambda[1] <- 1
  for (item in 2:I){
    lambda[item] ~ dnorm(lambda.mean.0, lambda.precision.0)
  }

  # prior distribution for the factor variance
  factor.precision ~ dgamma(factor.alpha.0, factor.beta.0)

  # saved parameters
  factor.variance <- 1/factor.precision
  for (item in 1:I){
    psi[item] <- 1/inv.psi[item]
    itemValue[item] <- lambda[item]*inv.psi[item]*lambda[item]
  }
  
  # reliability for sum scores:
  sumReliability <- sum(lambda[])*sum(lambda[]) / (sum(lambda[])*sum(lambda[]) + sum(psi[]))

  # reliability for factor scores:
  posteriorFactorVariance <- 1 / (factor.variance + sum(itemValue[]))
  factorReliability <- factor.variance / (factor.variance + posteriorFactorVariance)
}
"

# next, create data for JAGS to use:
model01.data = list(
  N = nrow(data),
  X = data,
  I = nItems,
  mu.mean.0 = mu.mean.0,
  mu.precision.0 = mu.precision.0,
  lambda.mean.0 = lambda.mean.0,
  lambda.precision.0 = lambda.precision.0,
  psi.alpha.0 = psi.alpha.0,
  psi.beta.0 = psi.beta.0,
  factor.alpha.0 = factor.alpha.0,
  factor.beta.0 = factor.beta.0
)


# for reproducable seeds (without parallel JAGS)
model01.seed = 27032019

#   Note: here, the random number seed cannot be the same per seed or the chains will be the same
RNGname = c("Wichmann-Hill","Marsaglia-Multicarry", "Super-Duper", "Mersenne-Twister")
if (RNGname[1] %in% c("Wichmann-Hill", "Marsaglia-Multicarry", 
                      "Super-Duper", "Mersenne-Twister")) {
  RNGname[1] <- paste("base::", RNGname[1], sep = "")
}

model01.init.values <- vector("list", nchains)
for (i in 1:nchains) {
  model01.init.values[[i]]$.RNG.name <- RNGname[1]
  model01.init.values[[i]]$.RNG.seed <- model01.seed + i
}


model01.parameters = c("mu", "lambda",  "psi", "factor.variance", "deviance", 
                       "posteriorFactorVariance", "factorReliability", "omega")

# load two JAGS modules: GLM (which makes the analysis go more efficiently) and DIC (which gives an index of relative fit)
load.module("glm")
load.module("dic")

# Submit model to JAGS for compiling
model01.JAGS = jags.model(
  file = textConnection(model01.syntax),
  data = model01.data,
  n.chains = nchains,
  n.adapt = nadapt,
  inits = model01.init.values
)

# Draw samples; NOTE: BURNIN IS INCUDED IN SAMPLES
model01.Samples = coda.samples(
  model = model01.JAGS,
  variable.names = model01.parameters,
  n.iter = niter,
  thin = nthin
)

# trim off the burnin of the chain
model01.Posterior = window(x = model01.Samples, start = nburnin + 1, end = niter)

# check convergence (here, we have a constant with the marker item so we need more syntax):
model01.gelmandiag = lapply(X = model01.Posterior, FUN = function(x) return(x[,which(colnames(x) != "lambda[1]")]))
gelman.diag(model01.gelmandiag)

# look at results 
summary(model01.Posterior)

# next: conduct ppmc

# list number of simulated data sets
nSimulatedDataSets = 5000

# create one large matrix of posterior value by disentangling chains
model01.Posterior.all = do.call("rbind", model01.Posterior)

# determine columns of posterior that go into each model matrix
muCols = grep(x = colnames(model01.Posterior.all), pattern = "mu")
lambdaCols = grep(x = colnames(model01.Posterior.all), pattern = "lambda")
psiCols = grep(x = colnames(model01.Posterior.all), pattern = "psi")
phiCols = grep(x = colnames(model01.Posterior.all), pattern = "factor.variance")

# save simulated correlations:
simCor = matrix(data = NA, nrow = nSimulatedDataSets, ncol = nItems*(nItems-1)/2)
simCovModel01b = matrix(data = NA, nrow = nSimulatedDataSets, ncol = nItems*nItems)
simSRMR = matrix(data = NA, nrow = nSimulatedDataSets, ncol = 1)
simRMSEA = matrix(data = NA, nrow = nSimulatedDataSets, ncol = 1)

# save model-based covariances:
model01Cov = matrix(data = NA, nrow = nSimulatedDataSets, ncol = nItems*nItems)

# model DF
modelDF = (nItems*(nItems+1)/2) - nItems*2

dataCov = cov(data)
detDataCov = det(dataCov)

# loop through data sets (can be sped up with functions and lapply)
pb = txtProgressBar()
sim = 1
for (sim in 1:nSimulatedDataSets){
  
  # draw sample from one iteration of posterior chain 
  iternum = sample(x = 1:nrow(model01.Posterior.all), size = 1, replace = TRUE)
  
  # get parameters for that sample: put into factor model matrices for easier generation of data
  mu = matrix(data = model01.Posterior.all[iternum, muCols], ncol = 1)
  lambda = matrix(data = model01.Posterior.all[iternum, lambdaCols], ncol = 1)
  psi = diag(model01.Posterior.all[iternum, psiCols])
  phi = matrix(model01.Posterior.all[iternum, phiCols])
  
  # create model-implied mean and covariance matrix (marginal for X)
  meanVec = mu
  covMat = lambda %*% phi %*% t(lambda) + psi
  model01Cov[sim, ] = c(covMat)
  
  # randomly draw data with same sample size from MVN with mean=meanVec and cov=covMat
  simData = rmvnorm(n = nrow(data), mean = meanVec, sigma = covMat)
  
  # create sample statistics from simulated data (we'll use correlation matrix, starting with upper triangle)
  simCor[sim,] = matrix(data = c(cor(simData)[upper.tri(cor(simData))]), nrow = 1)
  
  # calculate the value of SRMR using simulated data's covariance matrix and observed covariance matrix
  simCov = cov(simData)
  simCovModel01b[sim,] = c(cov(simData))
  difCov = dataCov-simCov
  stdDifCov = sqrt(solve(diag(diag(dataCov)))) %*% difCov %*% sqrt(solve(diag(diag(dataCov))))
  
  # using formula from book:
  simSRMR[sim,1] = sum(colSums(stdDifCov*stdDifCov))/((ncol(simCov)*(ncol(simCov)-1))/2)
  
  setTxtProgressBar(pb = pb, value = sim/nSimulatedDataSets)
}
close(pb)

# first, we examine the posterior predictive distribution of SRMR (p. 241)
hist(simSRMR[,1])
plot(density(simSRMR))
quantile(simSRMR)
mean(simSRMR)

# label values of simCor to ensure we have the right comparison
corNames = NULL
for (i in 1:(ncol(simData)-1)){
  for (j in (i+1):ncol(simData)){
    corNames = c(corNames, paste0("cor", i, "." , j))
  }
}
colnames(simCor) = corNames

# show how one correlation compares to distribution of simulated correlations
dataCor = cor(data)
hist(simCor[,1])
plot(density(simCor[,1]))
lines(x = c(dataCor[1,2], dataCor[1,2]), y = c(0, max(density(simCor[,1])$y)), col = 2)
quantile(simCor[,1])
mean(simCor[,1])

# create quantiles of correlations to see where each observed correlation falls
corQuantiles = NULL

# compute the quantiles of the observed correlations:

col = 1
for (i in 1:(ncol(simData)-1)){
  for (j in (i+1):ncol(simData)){
    # get empirical CDF of simulated correlation distribution
    corEcdf = ecdf(simCor[,col])
    corQuantiles = rbind(corQuantiles, c(i, j, summary(corEcdf), dataCor[i,j], corEcdf(dataCor[i,j])))
    
    col = col + 1
  }
}
colnames(corQuantiles)[1:2] = c("Item 1", "Item 2")
colnames(corQuantiles)[9:10] = c("ObsCor", "CorPctile")
corQuantiles[which(corQuantiles[,10] > .975 | corQuantiles[,10] < .025),]


# next: examine reliabilities
keepCols = which(colnames(model01.Posterior[[1]]) %in% c("sumReliability", "factorReliability"))
plot(mcmc.list(lapply(X = model01.Posterior, FUN = function(x) return(mcmc(x[,keepCols])))))
