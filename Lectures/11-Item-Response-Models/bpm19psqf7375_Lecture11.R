# Install/Load Packages ===============================================================================================
if (!require(R2jags)) install.packages("R2jags")
library(R2jags)

if (!require(CDM)) install.packages("CDM")
library(CDM)


# We will use the Tatsuoka (1984) fraction subtraction data for today's examples
# See DeCarlo (2011, p. 9) for the items: https://scholar.google.com/scholar?hl=en&as_sdt=0%2C36&q=l+decarlo+2011&btnG=



# First, we will treat these data as unidimensional to demonstrate unidimensional IRT models.
#   We can use the syntax from the unidimensional CFA model as a start for modeling the 2PL model:
#     This uses slope/intercept form, which we will change to discrimination/difficulty later

# model specs:
nItems = ncol(FSdata)
nchains = 4
niter = 10000
nburnin = 2000
nadapt = 2000
nthin = 1


# marker item:
model01.syntax = "
model{
  # measurement model specification
    for (person in 1:N){
      for (item in 1:I){
        X[person, item] ~ dbern(phi(mu[item] + lambda[item]*xfactor[person]))
      }
    }

  # prior distributions for the factor:
    for (person in 1:N){
      xfactor[person] ~ dnorm(0, factor.precision)
    }

  # prior distributions for the measurement model mean/precision parameters
    for (item in 1:I){
      mu[item] ~ dnorm(mu.mean.0, mu.precision.0)
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
}


"

# specification of prior values for measurement model parameters:
#   item intercepts
mu.mean.0 = 0
mu.variance.0 = 1000000
mu.precision.0 = 1 / mu.variance.0

#   Factor loadings -- these are the discriminations
lambda.mean.0 = 0
lambda.variance.0 = 1000000
lambda.precision.0 = 1 / lambda.variance.0

# unique variances -- these do not exist

# values for prior for factor variance (based on variance of marker item)
factor.df.0 = 1
factor.var.0 = 1
factor.alpha.0 = factor.df.0/2
factor.beta.0 = (factor.df.0*factor.var.0)/2


# next, create data for JAGS to use:
model01.data = list(
  N = nrow(FSdata),
  X = FSdata,
  I = nItems,
  mu.mean.0 = mu.mean.0,
  mu.precision.0 = mu.precision.0,
  lambda.mean.0 = lambda.mean.0,
  lambda.precision.0 = lambda.precision.0,
  factor.alpha.0 = factor.alpha.0,
  factor.beta.0 = factor.beta.0
)

# copy initial values

# for reproducable seeds (without parallel JAGS)
model01.seed = 06042019

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

# for xfactor, we are only going to use
model01.parameters = c("mu", "lambda",  "factor.variance", "deviance", "xfactor[23]", "xfactor[28]", "xfactor[55]")

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

# Look at sampling algorithms
list.samplers(model01.JAGS)

# Draw samples; NOTE: BURNIN IS INCUDED IN SAMPLES
model01.JAGS.Samples = coda.samples(
  model = model01.JAGS,
  variable.names = model01.parameters,
  n.iter = niter,
  thin = nthin
)

save.image()
# trim off the burnin of the chain
model01.Posterior = window(x = model01.JAGS.Samples, start = nburnin + 1, end = niter)

# check convergence (here, we have a constant with the marker item so we need more syntax):
model01.gelmanDiag = lapply(X = model01.Posterior, FUN = function(x) return(x[,which(colnames(x) != "lambda[1]")]))
gelman.diag(model01.gelmanDiag)

# look at results 
plot(mcmc.list(lapply(X = model01.Posterior, FUN = function(x) return(x[,which(colnames(x) == "factor.variance")]))))
summary(model01.Posterior)


# lets compare deviance from model01b to model01a:
summary(model01a.Posterior)




# marker item:
model01.syntax = "
model{
# measurement model specification
for (person in 1:N){
for (item in 1:I){
prob[person, item] = 1/(1+exp(-1*(mu[item] + lambda[item]*xfactor[person])))
X[person, item] ~ dbern(prob[person,item])    
}
}

# prior distributions for the factor:
for (person in 1:N){
xfactor[person] ~ dnorm(0, factor.precision)
}

# prior distributions for the measurement model mean/precision parameters
for (item in 1:I){
mu[item] ~ dnorm(mu.mean.0, mu.precision.0)
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
}


"

# next, create data for JAGS to use:
model01.data = list(
  N = nrow(FSdata),
  X = FSdata,
  I = nItems,
  mu.mean.0 = mu.mean.0,
  mu.precision.0 = mu.precision.0,
  lambda.mean.0 = lambda.mean.0,
  lambda.precision.0 = lambda.precision.0,
  factor.alpha.0 = factor.alpha.0,
  factor.beta.0 = factor.beta.0
)

# copy initial values

# for reproducable seeds (without parallel JAGS)
model01.seed = 06042019

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

#, "xfactor") --- note, we'll leave out the factor scores for now. These will come back once we know the model works.
model01.parameters = c("mu", "lambda",  "factor.variance", "deviance", "xfactor")

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

# Look at sampling algorithms
list.samplers(model01.JAGS)

# Draw samples; NOTE: BURNIN IS INCUDED IN SAMPLES
model01.JAGS.Samples = coda.samples(
  model = model01.JAGS,
  variable.names = model01.parameters,
  n.iter = niter,
  thin = nthin
)

# trim off the burnin of the chain
model01.Posterior = window(x = model01.JAGS.Samples, start = nburnin + 1, end = niter)

# check convergence (here, we have a constant with the marker item so we need more syntax):
model01.gelmanDiag = lapply(X = model01.Posterior, FUN = function(x) return(x[,which(colnames(x) != "lambda[1]")]))
gelman.diag(model01.gelmanDiag)

# look at results 
plot(mcmc.list(lapply(X = model01.Posterior, FUN = function(x) return(x[,which(colnames(x) == "factor.variance")]))))
summary(model01.Posterior)



