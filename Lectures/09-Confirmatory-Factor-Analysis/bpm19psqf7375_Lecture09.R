

# Install/Load Packages ===============================================================================================
if (!require(rjags)) install.packages("rjags")
library(rjags)

if (!require(mvtnorm)) install.packages("mvtnorm")
library(mvtnorm)

# Import/Examine Data =================================================================================================
# Today's example is from a bootstrap resample of 177 undergradutes at a large state university in the midwest. The 
# survey was a measure of 10 questions about their beliefs in various conspiracy theories that were being passed around 
# the internet in the early 2010s. Additionally, gender was included in the survey. All items responses were on a 5-
# point Likert scale where 1=Strongly Disagree, 2=Disagree, 3=Neither Agree or Disagree, 4=Agree, and 5=Strongly Agree.

# Questions:
# 1. The U.S. invasion of Iraq was not part of a campaign to fight terrorism, but was driven by oil companies and Jews 
#     in the U.S. and Israel.
# 2. Certain U.S. government officials planned the attacks of September 11, 2001 because they wanted the United States
#     to go to war in the Middle East.
# 3. President Barack Obama was not really born in the United States and does not have an authentic Hawaiian birth 
#     certificate.
# 4. The current financial crisis was secretly orchestrated by a small group of Wall Street bankers to extend the power 
#     of the Federal Reserve and further their control of the world's economy.
# 5. Vapor trails left by aircraft are actually chemical agents deliberately sprayed in a clandestine program directed 
#     by government officials.
# 6. Billionaire George Soros is behind a hidden plot to destabilize the American government, take control of the 
#     media, and put the world under his control.
# 7. The U.S. government is mandating the switch to compact fluorescent light bulbs because such lights make people 
#     more obedient and easier to control.
# 8. Government officials are covertly Building a 12-lane \"NAFTA superhighway\" that runs from Mexico to Canada 
#     through America's heartland.
# 9. Government officials purposely developed and spread drugs like crack-cocaine and diseases like AIDS in order to 
#     destroy the African American community.
# 10. God sent Hurricane Katrina to punish America for its sins.

# read in data:
conspiracy = read.csv("Lectures/09-Confirmatory-Factor-Analysis/conspiracies.csv")

# plot each item
hist(conspiracy$PolConsp1, main = "PolConsp1", xlab = "1. The U.S. invasion of Iraq was not part of a campaign to fight terrorism, but was driven by oil companies and Jews in the U.S. and Israel.")
hist(conspiracy$PolConsp2, main = "PolConsp2", xlab = "2. Certain U.S. government officials planned the attacks of September 11, 2001 because they wanted the United States to go to war in the Middle East.")
hist(conspiracy$PolConsp3, main = "PolConsp3", xlab = "3. President Barack Obama was not really born in the United States and does not have an authentic Hawaiian birth certificate.")
hist(conspiracy$PolConsp4, main = "PolConsp4", xlab = "4. The current financial crisis was secretly orchestrated by a small group of Wall Street bankers to extend the power of the Federal Reserve and further their control of the world's economy.")
hist(conspiracy$PolConsp5, main = "PolConsp5", xlab = "5. Vapor trails left by aircraft are actually chemical agents deliberately sprayed in a clandestine program directed by government officials.")
hist(conspiracy$PolConsp6, main = "PolConsp6", xlab = "6. Billionaire George Soros is behind a hidden plot to destabilize the American government, take control of the media, and put the world under his control.")
hist(conspiracy$PolConsp7, main = "PolConsp7", xlab = "7. The U.S. government is mandating the switch to compact fluorescent light bulbs because such lights make people more obedient and easier to control.")
hist(conspiracy$PolConsp8, main = "PolConsp8", xlab = "8. Government officials are covertly Building a 12-lane \"NAFTA superhighway\" that runs from Mexico to Canada through America's heartland.")
hist(conspiracy$PolConsp9, main = "PolConsp9", xlab = "9. Government officials purposely developed and spread drugs like crack-cocaine and diseases like AIDS in order to destroy the African American community.")
hist(conspiracy$PolConsp10, main = "PolConsp10", xlab = "10. God sent Hurricane Katrina to punish America for its sins.")

# model 1a: single factor model with standardized factor scale identification ==========================================

# model specs:
nItems = ncol(conspiracy[paste0("PolConsp", 1:10)])
nchains = 4
niter = 3000
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
psi.var.0 = apply(X = conspiracy[paste0("PolConsp", 1:10)], MARGIN = 2, FUN = var)
psi.alpha.0 = psi.df.0 / 2
psi.beta.0 = (psi.df.0 * psi.var.0) / 2

model01a.data = list(
  N = nrow(conspiracy),
  X = conspiracy[paste0("PolConsp", 1:10)],
  I = nItems,
  mu.mean.0 = mu.mean.0,
  mu.precision.0 = mu.precision.0,
  lambda.mean.0 = lambda.mean.0,
  lambda.precision.0 = lambda.precision.0,
  psi.alpha.0 = psi.alpha.0,
  psi.beta.0 = psi.beta.0
)

#, "xfactor") --- note, we'll leave out the factor scores for now. These will come back once we know the model works.
model01a.parameters = c("mu", "lambda",  "psi", "deviance")


# for reproducable seeds (without parallel JAGS)
model01a.seed = 23022019

#   Note: here, the random number seed cannot be the same per seed or the chains will be the same
RNGname = c("Wichmann-Hill","Marsaglia-Multicarry", "Super-Duper", "Mersenne-Twister")
if (RNGname[1] %in% c("Wichmann-Hill", "Marsaglia-Multicarry", 
                      "Super-Duper", "Mersenne-Twister")) {
  RNGname[1] <- paste("base::", RNGname[1], sep = "")
}

model01a.init.values <- vector("list", nchains)
for (i in 1:nchains) {
  model01a.init.values[[i]]$.RNG.name <- RNGname[1]
  model01a.init.values[[i]]$.RNG.seed <- model01a.seed + i
}


model01a.syntax = "
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
  xfactor[person] ~ dnorm(0, 1)
}

# prior distributions for the measurement model parameters
for (item in 1:I){
  mu[item] ~ dnorm(mu.mean.0, mu.precision.0)
  lambda[item] ~ dnorm(lambda.mean.0, lambda.precision.0)
  inv.psi[item] ~ dgamma(psi.alpha.0, psi.beta.0[item])
}

# saved parameters
for (item in 1:I){
  psi[item] <- 1/inv.psi[item]
}
}
"


# load two JAGS modules: GLM (which makes the analysis go more efficiently) and DIC (which gives an index of relative fit)
load.module("glm")
load.module("dic")

# Submit model to JAGS for compiling. We'll turn adaptation off to control it in the next line of syntax:
model01a.JAGS = jags.model(
  file = textConnection(model01a.syntax),
  data = model01a.data,
  n.chains = nchains,
  n.adapt = nadapt,
  inits = model01a.init.values
)

list.samplers(model01a.JAGS)

# Draw samples; NOTE: BURNIN IS INCUDED IN SAMPLES
model01a.Samples = coda.samples(
  model = model01a.JAGS,
  variable.names = model01a.parameters,
  n.iter = niter,
  thin = nthin
)



gelman.diag(model01a.Samples)
summary(model01a.Samples)

# the factor loadings are not in good shape...
# plot(model01a.Samples)

# even after we trim off the burnin period:
model01a.Posterior = window(x = model01a.Samples, start = nburnin + 1, end = niter)
gelman.diag(model01a.Posterior)

# plot(model01a.Posterior)

# Model 1b: Single Factor Model with Marker-Item Factor Variance Scale Identification =================================

# model specs:
nItems = ncol(conspiracy[paste0("PolConsp", 1:10)])
nchains = 4
niter = 3000
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
psi.var.0 = apply(X = conspiracy[paste0("PolConsp", 1:10)], MARGIN = 2, FUN = var)
psi.alpha.0 = psi.df.0 / 2
psi.beta.0 = (psi.df.0 * psi.var.0) / 2

# values for prior for factor variance (based on variance of marker item)
factor.df.0 = 1
factor.var.0 = var(conspiracy$PolConsp1)
factor.alpha.0 = factor.df.0/2
factor.beta.0 = (factor.df.0*factor.var.0)/2

# marker item:
model01b.marker.syntax = "
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
  inv.psi[item] ~ dgamma(psi.alpha.0, psi.beta.0[item])
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
}
}
"

# next, create data for JAGS to use:
model01b.data = list(
  N = nrow(conspiracy),
  X = conspiracy[paste0("PolConsp", 1:10)],
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

# copy initial values
model01b.init.values = model01a.init.values

#, "xfactor") --- note, we'll leave out the factor scores for now. These will come back once we know the model works.
model01b.parameters = c("mu", "lambda",  "psi", "factor.variance", "deviance")

# load two JAGS modules: GLM (which makes the analysis go more efficiently) and DIC (which gives an index of relative fit)
load.module("glm")
load.module("dic")

# Submit model to JAGS for compiling
model01b.marker.JAGS = jags.model(
  file = textConnection(model01b.marker.syntax),
  data = model01b.data,
  n.chains = nchains,
  n.adapt = nadapt,
  inits = model01b.init.values
)

list.samplers(model01b.marker.JAGS)

# Draw samples; NOTE: BURNIN IS INCUDED IN SAMPLES
model01b.marker.Samples = coda.samples(
  model = model01b.marker.JAGS,
  variable.names = model01b.parameters,
  n.iter = niter,
  thin = nthin
)

# trim off the burnin of the chain
model01b.marker.Posterior = window(x = model01b.marker.Samples, start = nburnin + 1, end = niter)

# convergence check fails: we have a constant with the marker item
gelman.diag(model01b.marker.Posterior)

# check convergence (here, we have a constant with the marker item so we need more syntax):
model01b.gelmandiag = lapply(X = model01b.marker.Posterior, FUN = function(x) return(x[,which(colnames(x) != "lambda[1]")]))
gelman.diag(model01b.gelmandiag)

# look at results 
summary(model01b.marker.Posterior)

# now things look converged!

# lets compare deviance from model01b to model01a:
summary(model01a.Posterior)

# almost the same...caution: deviance can be misleading for overall chain convergence

# model is converged, but does the model fit the data?

# Posterior Predictive Checks for Model 1b ============================================================================

# list number of simulated data sets
nSimulatedDataSets = 5000

# create one large matrix of posterior value by disentangling chains
model01b.Posterior.all = do.call("rbind", model01b.marker.Posterior)

# determine columns of posterior that go into each model matrix
muCols = grep(x = colnames(model01b.Posterior.all), pattern = "mu")
lambdaCols = grep(x = colnames(model01b.Posterior.all), pattern = "lambda")
psiCols = grep(x = colnames(model01b.Posterior.all), pattern = "psi")
phiCols = grep(x = colnames(model01b.Posterior.all), pattern = "factor.variance")

# save simulated correlations:
simCor = matrix(data = NA, nrow = nSimulatedDataSets, ncol = nItems*(nItems-1)/2)
simCovModel01b = matrix(data = NA, nrow = nSimulatedDataSets, ncol = nItems*nItems)
simSRMR = matrix(data = NA, nrow = nSimulatedDataSets, ncol = 1)
simRMSEA = matrix(data = NA, nrow = nSimulatedDataSets, ncol = 1)

# save model-based covariances:
model01bCov = matrix(data = NA, nrow = nSimulatedDataSets, ncol = nItems*nItems)

# model DF
modelDF = (nItems*(nItems+1)/2) - 20

dataCov = cov(conspiracy[paste0("PolConsp", 1:10)])
detDataCov = det(dataCov)

# loop through data sets (can be sped up with functions and lapply)
pb = txtProgressBar()
sim = 1
for (sim in 1:nSimulatedDataSets){
  
  # draw sample from one iteration of posterior chain 
  iternum = sample(x = 1:nrow(model01b.Posterior.all), size = 1, replace = TRUE)
  
  # get parameters for that sample: put into factor model matrices for easier generation of data
  mu = matrix(data = model01b.Posterior.all[iternum, muCols], ncol = 1)
  lambda = matrix(data = model01b.Posterior.all[iternum, lambdaCols], ncol = 1)
  psi = diag(model01b.Posterior.all[iternum, psiCols])
  phi = matrix(model01b.Posterior.all[iternum, phiCols])
  
  # create model-implied mean and covariance matrix (marginal for X)
  meanVec = mu
  covMat = lambda %*% phi %*% t(lambda) + psi
  model01bCov[sim, ] = c(covMat)
    
  # randomly draw data with same sample size from MVN with mean=meanVec and cov=covMat
  simData = rmvnorm(n = nrow(conspiracy), mean = meanVec, sigma = covMat)
  
  # create sample statistics from simulated data (we'll use correlation matrix, starting with upper triangle)
  simCor[sim,] = matrix(data = c(cor(simData)[upper.tri(cor(simData))]), nrow = 1)
  
  # calculate the value of SRMR using simulated data's covariance matrix and observed covariance matrix
  simCov = cov(simData)
  simCovModel01b[sim,] = c(cov(simData))
  difCov = dataCov-simCov
  stdDifCov = sqrt(solve(diag(diag(dataCov)))) %*% difCov %*% sqrt(solve(diag(diag(dataCov))))

  # using formula from book:
  simSRMR[sim,1] = sum(colSums(stdDifCov*stdDifCov))/((ncol(simCov)*(ncol(simCov)-1))/2)
  
  # can also do a similar process for RMSEA (assuming covariance matrix is from ML estimation - discrepancy function)
  simRMSEA[sim,1] = sqrt((log(det(simCov)) - log(det(dataCov)) + sum(diag(dataCov %*% solve(simCov))) - nItems)/modelDF)
  
  setTxtProgressBar(pb = pb, value = sim/nSimulatedDataSets)
}
close(pb)

# first, we examine the posterior predictive distribution of SRMR (p. 241)
hist(simSRMR[,1])
plot(density(simSRMR))
quantile(simSRMR)
mean(simSRMR)

# next we can examine the posterior predictive distribution of RMSEA
hist(simRMSEA[,1])
plot(density(simRMSEA))
quantile(simRMSEA)
mean(simRMSEA)

# label values of simCor to ensure we have the right comparison
corNames = NULL
for (i in 1:(ncol(simData)-1)){
  for (j in (i+1):ncol(simData)){
    corNames = c(corNames, paste0("cor", i, "." , j))
  }
}
colnames(simCor) = corNames

# show how one correlation compares to distribution of simulated correlations
dataCor = cor(conspiracy[paste0("PolConsp", 1:10)])
hist(simCor[,1])
plot(density(simCor[,1]))
lines(x = c(dataCor[1,2], dataCor[1,2]), y = c(0, 5), col = 2)
quantile(simRMSEA[,1])
mean(simRMSEA)

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

# Model 2a: Saturated Mean/Variance/Covariance Model ==================================================================

# model specs:
nItems = ncol(conspiracy[paste0("PolConsp", 1:10)])
nchains = 4
niter = 3000
nburnin = 2000
nadapt = 2000
nthin = 1


# prior values for each mean -- joinly MVN (but diagonal covariance/precision matrix)
mean.0 = 0
mean.variance.0 = 1000000
mean.precision.0 = 1/mean.variance.0
meanVec.0 = matrix(data = rep(mean.0, nItems), nrow = nItems, ncol = 1)
meanPrecision.0 = diag(nItems)*mean.precision.0

# prior values for variances -- variances of observed data
R.0 = (apply(X = conspiracy[paste0("PolConsp", 1:10)], MARGIN = 2, FUN = var) * diag(nItems))/nItems
k.0 = nItems

# load data into list
model02a.data = list(N = nrow(conspiracy),
                     X = conspiracy[paste0("PolConsp", 1:10)],
                     I = nItems, 
                     meanVec.0 = meanVec.0,
                     meanPrecision.0 = meanPrecision.0,
                     R.0 = R.0,
                     k.0 = k.0)

# save parameters
model02a.parameters = c("mu", "sigma", "deviance")

model02a.syntax = "
model{
  # define model likelihood as multivariate normal
  for (person in 1:N){
    X[person, 1:I] ~ dmnorm(mu[1:I], sigmainv[1:I, 1:I])
  }

  # prior for mean vector
  mu[1:I] ~ dmnorm(meanVec.0, meanPrecision.0) 

  # prior for inverse covariance matrix
  sigmainv[1:I, 1:I] ~ dwish(R.0, k.0)

  # save covariance matrix 
  sigma = inverse(sigmainv)

}
"


# copy initial values
model02a.init.values = model01a.init.values

# load two JAGS modules: GLM (which makes the analysis go more efficiently) and DIC (which gives an index of relative fit)
load.module("glm")
load.module("dic")

# Submit model to JAGS for compiling. We'll turn adaptation off to control it in the next line of syntax:
model02a.JAGS = jags.model(
  file = textConnection(model02a.syntax),
  data = model02a.data,
  n.chains = nchains,
  n.adapt = nadapt,
  inits = model02a.init.values
)



list.samplers(model02a.JAGS)

# Draw samples; NOTE: BURNIN IS INCUDED IN SAMPLES
model02a.Samples = coda.samples(
  model = model02a.JAGS,
  variable.names = model02a.parameters,
  n.iter = niter,
  thin = nthin
)

# trim off the burnin of the chain
model02a.Posterior = window(x = model02a.Samples, start = nburnin + 1, end = niter)

# convergence check fails: there is an issue with the multivariate version
gelman.diag(model02a.Posterior)
gelman.diag(model02a.Posterior, multivariate = FALSE)
# plot(model02a.Posterior)

summary(model02a.Posterior)

# posterior predictive check of model 2a =============================================================================================

# let's look at the distribution of the covariances from this model, as compared to the (ML) data covariances and simulated ones

# compute the quantiles of the chain covariances:
dataCov = cov(conspiracy[paste0("PolConsp", 1:10)])
chainCovs = mcmc.list(lapply(X = model02a.Posterior, FUN = function(x) return(x[,12:ncol(x)])))
chainCovsS = summary(chainCovs)
chainCovQ = cbind(chainCovsS$quantiles, chainCovsS$statistics, matrix(data = NA, nrow = nrow(chainCovsS$statistics), ncol = 2))
row = 1
for (i in 1:ncol(dataCov)){
  for (j in 1:ncol(dataCov)){
    chainCovQ[row, 10] = dataCov[i,j]
    chainCovQ[row, 11] = chainCovQ[row, 6] - chainCovQ[row, 10]
    row = row +1
  }
}
colnames(chainCovQ)[10:11] = c("dataCov", "DifDataCovEAP")
chainCovQ

# comparison of model vs. chain covariance for all pairs of items
plot(chainCovQ[, 6], chainCovQ[, 10], ylab = "dataCov", xlab = "modelCov")

# now, let's compare this with the model-implied covariance matrix from model01b
model02a.Posterior.all = do.call("rbind", chainCovs)
plot(density(model02a.Posterior.all[,1]), ylim = c(0, max(density(model01bCov[,1])$y, density(model02a.Posterior.all[,1])$y)+ .1),
     main = "Item 1 Variance")
lines(density(model01bCov[,1]), col = 2, lty = 2)

plot(density(model02a.Posterior.all[,2]), ylim = c(0, max(density(model01bCov[,2])$y, density(model02a.Posterior.all[,2])$y)+ .1),
     main = "Item 1 and Item 2 Covariance")
lines(density(model01bCov[,2]), col = 2, lty = 2)

# we can also generate a posterior predictive distribution from this model as well:

# list number of simulated data sets
nSimulatedDataSets = 5000

# create one large matrix of posterior value by disentangling chains
model02a.Posterior.all2 = do.call("rbind", model02a.Posterior)

# determine columns of posterior that go into each model matrix
muCols = grep(x = colnames(model02a.Posterior.all2), pattern = "mu")
sigmaCols = grep(x = colnames(model02a.Posterior.all2), pattern = "sigma")

# save simulated correlations:
simCorModel02a = matrix(data = NA, nrow = nSimulatedDataSets, ncol = nItems*(nItems-1)/2)
simCovModel02a = matrix(data = NA, nrow = nSimulatedDataSets, ncol = nItems*nItems)

# loop through data sets (can be sped up with functions and lapply)
pb = txtProgressBar()
sim = 1
for (sim in 1:nSimulatedDataSets){
  
  # draw sample from one iteration of posterior chain 
  iternum = sample(x = 1:nrow(model02a.Posterior.all2), size = 1, replace = TRUE)
  
  # get parameters for that sample: put into factor model matrices for easier generation of data
  mu = matrix(data = model02a.Posterior.all2[iternum, muCols], ncol = 1)
  sigma = matrix(data = model02a.Posterior.all2[iternum, sigmaCols], ncol = 10, nrow = 10)
  
  # create model-implied mean and covariance matrix (marginal for X)
  meanVec = mu
  covMat = sigma
  
  # randomly draw data with same sample size from MVN with mean=meanVec and cov=covMat
  simData = rmvnorm(n = nrow(conspiracy), mean = meanVec, sigma = covMat)
  
  # create sample statistics from simulated data (we'll use correlation matrix, starting with upper triangle)
  simCorModel02a[sim,] = matrix(data = c(cor(simData)[upper.tri(cor(simData))]), nrow = 1)
  
  # calculate the value of SRMR using simulated data's covariance matrix and observed covariance matrix
  simCovModel02a[sim,] = c(cov(simData))
  # difCov = dataCov-simCov
  # stdDifCov = sqrt(solve(diag(diag(dataCov)))) %*% difCov %*% sqrt(solve(diag(diag(dataCov))))
  
  # using formula from book:
  # simSRMR[sim,1] = sum(colSums(stdDifCov*stdDifCov))/((ncol(simCov)*(ncol(simCov)-1))/2)
  
  # can also do a similar process for RMSEA (assuming covariance matrix is from ML estimation - discrepancy function)
  # simRMSEA[sim,1] = sqrt((log(det(simCov)) - log(det(dataCov)) + sum(diag(dataCov %*% solve(simCov))) - nItems)/modelDF)
  
  setTxtProgressBar(pb = pb, value = sim/nSimulatedDataSets)
}

# next, let's see how the simulated covariances compare with the model implied ones (from Model02a only)

# now, let's compare this with the model-implied covariance matrix from model01b
plot(density(model02a.Posterior.all[,1]), ylim = c(0, max(density(simCovModel02a[,1])$y, density(model02a.Posterior.all[,1])$y) + .1),
     main = "Item 1 Variance")
lines(density(simCovModel02a[,1]), col = 2, lty = 2)

plot(density(model02a.Posterior.all[,2]), ylim = c(0, max(density(simCovModel02a[,2])$y, density(model02a.Posterior.all[,2])$y)+ .1),
     main = "Item 1 and Item 2 Covariance")
lines(density(simCovModel02a[,2]), col = 2, lty = 2)

# now comparing the simulated correlations between models 1b and 2a
plot(density(simCovModel01b[,1]), ylim = c(0, max(density(simCovModel02a[,1])$y, density(simCovModel01b[,1])$y) + .1),
     main = "Item 1 Variance")
lines(density(simCovModel02a[,1]), col = 2, lty = 2)

qqplot(y = simCovModel01b[,1], x = simCovModel02a[,1], main = "Q-Q plot of Item 1 and Item 2 Covariance", ylab = "Model01b", xlab = "model02a",
       xlim = c(min(c(simCovModel01b[,1], simCovModel02a[,1])), max(c(simCovModel01b[,1], simCovModel02a[,1]))),
       ylim = c(min(c(simCovModel01b[,1], simCovModel02a[,1])), max(c(simCovModel01b[,1], simCovModel02a[,1]))))
lines(x = c(min(c(simCovModel01b[,1], simCovModel02a[,1])), max(c(simCovModel01b[,1], simCovModel02a[,1]))), 
      y = c(min(c(simCovModel01b[,1], simCovModel02a[,1])), max(c(simCovModel01b[,1], simCovModel02a[,1]))))
ks.test(simCovModel01b[,1], simCovModel02a[,1])


plot(density(simCovModel01b[,2]), ylim = c(0, max(density(simCovModel02a[,2])$y, density(simCovModel01b[,2])$y)+ .1),
     main = "Item 1 and Item 2 Covariance")
lines(density(simCovModel02a[,2]), col = 2, lty = 2)
qqplot(y = simCovModel01b[,2], x = simCovModel02a[,2], main = "Q-Q plot of Item 1 and Item 2 Covariance", ylab = "Model01b", xlab = "model02a",
       xlim = c(min(c(simCovModel01b[,2], simCovModel02a[,2])), max(c(simCovModel01b[,2], simCovModel02a[,2]))),
       ylim = c(min(c(simCovModel01b[,2], simCovModel02a[,2])), max(c(simCovModel01b[,2], simCovModel02a[,2]))))
lines(x = c(min(c(simCovModel01b[,2], simCovModel02a[,2])), max(c(simCovModel01b[,2], simCovModel02a[,2]))), 
      y = c(min(c(simCovModel01b[,2], simCovModel02a[,2])), max(c(simCovModel01b[,2], simCovModel02a[,2]))))
ks.test(simCovModel01b[,2], simCovModel02a[,2])


# Model Comparison -- can use DIC.samples to get DIC for both models ==================================================
#     Warning: The dic.samples function of rjags takes an extremely long time to run for the saturated model. 
#               So, we are switching back to R2jags, although the results are not repoducable exactly.

# getting DIC for our one-factor model:
model01b.function = function(){
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
    inv.psi[item] ~ dgamma(psi.alpha.0, psi.beta.0[item])
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
  }
}


model01b.r2jags = jags(
  data = model01b.data,
  inits = model01b.init.values,
  parameters.to.save = model01b.parameters,
  model.file = model01b.function,
  n.chains = nchains,
  n.iter = niter,
  n.burnin = nburnin,
  n.thin = nthin
)
model01b.r2jags

model02.function = function(){
  # define model likelihood as multivariate normal
  for (person in 1:N){
    X[person, 1:I] ~ dmnorm(mu[1:I], sigmainv[1:I, 1:I])
  }
  
  # prior for mean vector
  mu[1:I] ~ dmnorm(meanVec.0, meanPrecision.0) 
  
  # prior for inverse covariance matrix
  sigmainv[1:I, 1:I] ~ dwish(R.0, k.0)
  
  # save covariance matrix 
  sigma = inverse(sigmainv)
}


model02a.r2jags =  jags(
  data = model02a.data,
  inits = model02a.init.values,
  parameters.to.save = model02a.parameters,
  model.file = model02.function,
  n.chains = nchains,
  n.iter = niter,
  n.burnin = nburnin,
  n.thin = nthin
)
model02a.r2jags


# This syntax runs the example using dic.samples (but takes an excessive amount of time) 

# model01b.DICSamples = dic.samples(
#   model = model01b.marker.JAGS,
#   variable.names = model01b.parameters,
#   n.iter = niter,
#   thin = nthin,
#   type = "pD"
# )
# model01b.DICSamples

# 
# model02a.DICSamples = dic.samples(
#   model = model02a.JAGS,
#   variable.names = model02a.parameters,
#   n.iter = niter,
#   thin = nthin,
#   type = "pD"
# )
# model02a.DICSamples
# 
# diffdic(model01b.DICSamples, model02a.DICSamples)

# Model 03: A Two-factor model ========================================================================================

# model specs:
nItems = ncol(conspiracy[paste0("PolConsp", 1:10)])
nchains = 4
niter = 3000
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
psi.var.0 = apply(X = conspiracy[paste0("PolConsp", 1:10)], MARGIN = 2, FUN = var)
psi.alpha.0 = psi.df.0 / 2
psi.beta.0 = (psi.df.0 * psi.var.0) / 2

# values for prior for factor variance (based on variance of marker item)
factor.cov.0 = diag(2)
# factor.cov.0[1,2] = factor.cov.0[2,1] = .6
factor.invcov.0 = solve(factor.cov.0)
factor.invcov.df.0 = 2

# factor 1: general conspiracy: 1 3 4 6 10
# factor 2: government conspiracy: 2 5 7 8 9
model03a.function = function(){
  
  # measurement model specification
  for (person in 1:N){
      mean[person,  1] <- mu[1]  + lambda[1,  1]*xfactor[person, 1]
      mean[person,  2] <- mu[2]  + lambda[2,  2]*xfactor[person, 2]
      mean[person,  3] <- mu[3]  + lambda[3,  1]*xfactor[person, 1]
      mean[person,  4] <- mu[4]  + lambda[4,  1]*xfactor[person, 1]
      mean[person,  5] <- mu[5]  + lambda[5,  2]*xfactor[person, 2]
      mean[person,  6] <- mu[6]  + lambda[6,  1]*xfactor[person, 1]
      mean[person,  7] <- mu[7]  + lambda[7,  2]*xfactor[person, 2]
      mean[person,  8] <- mu[8]  + lambda[8,  2]*xfactor[person, 2]
      mean[person,  9] <- mu[9]  + lambda[9,  2]*xfactor[person, 2]
      mean[person, 10] <- mu[10] + lambda[10, 1]*xfactor[person, 1]
      
      for (item in 1:I){
        X[person, item] ~ dnorm(mean[person,item], inv.psi[item])  
      }
  }
  
  # prior distributions for the factor:
  for (person in 1:N){
    xfactor[person, 1:2] ~ dmnorm(kappa, inv.phi[,])
  }
  
  # fix factor means
  for (factor in 1:2){
    kappa[factor] <- 0
  }
  
  # prior distributions for the measurement model mean/precision parameters
  for (item in 1:I){
    mu[item] ~ dnorm(mu.mean.0, mu.precision.0)
    inv.psi[item] ~ dgamma(psi.alpha.0, psi.beta.0[item])
  }
  
  # prior distributions for the loadings (except the first loading, which is fixed to 1.0)
  lambda[1,1] <- 1
  lambda[3,1] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[4,1] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[6,1] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[10,1] ~ dnorm(lambda.mean.0, lambda.precision.0)
  lambda[2,2] <- 1
  lambda[5,2] ~ dnorm(lambda.mean.0, .1)
  lambda[7,2] ~ dnorm(lambda.mean.0, .1)
  lambda[8,2] ~ dnorm(lambda.mean.0, .1)
  lambda[9,2] ~ dnorm(lambda.mean.0, .1)
  
  # prior distribution for the factor covariance matrix
  inv.phi ~ dwish(factor.invcov.0, factor.invcov.df.0)
  
  # saved parameters

  for (item in 1:I){
    psi[item] <- 1/inv.psi[item]
  }
  phi = inverse(inv.phi)
  
}


# next, create data for JAGS to use:
model03a.data = list(
  N = nrow(conspiracy),
  X = conspiracy[paste0("PolConsp", 1:10)],
  I = nItems,
  mu.mean.0 = mu.mean.0,
  mu.precision.0 = mu.precision.0,
  lambda.mean.0 = lambda.mean.0,
  lambda.precision.0 = lambda.precision.0,
  psi.alpha.0 = psi.alpha.0,
  psi.beta.0 = psi.beta.0,
  factor.invcov.0 = factor.invcov.0,
  factor.invcov.df.0 = factor.invcov.df.0
)

# copy initial values
model03a.init.values = model01a.init.values

#, "xfactor") --- note, we'll leave out the factor scores for now. These will come back once we know the model works.
model03a.parameters = c("mu", "lambda",  "psi", "phi", "xfactor")

model03a.r2jags =  jags(
  data = model03a.data,
  inits = model03a.init.values,
  parameters.to.save = model03a.parameters,
  model.file = model03a.function,
  n.chains = nchains,
  n.iter = niter,
  n.burnin = nburnin,
  n.thin = nthin
)
model03a.r2jags

# examining the factor scores 

# Making It Easier ?: blavaan package ===================================================================================
if (!require(blavaan)) install.packages("blavaan")
library(blavaan)

# running the one factor model:
onefactor.blavaan.syntax = "
  factor =~ PolConsp1 + PolConsp2 + PolConsp3 + PolConsp4 + PolConsp5 +
            PolConsp6 + PolConsp7 + PolConsp8 + PolConsp9 + PolConsp10
  factor ~~ factor
"

onefactor.blavaan.analysis = bcfa(onefactor.blavaan.syntax, data = conspiracy)
summary(onefactor.blavaan.analysis) 

