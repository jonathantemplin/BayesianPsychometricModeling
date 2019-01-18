
# to run this script, you will first have to install JAGS from http://mcmc-jags.sourceforge.net

dataIQ = read.csv(file = "iqdata.csv")

#AUTOMATING PACKAGES NEEDED FOR ANALYSES--------------------------------------------------------------------
haspackage = require("R2jags")
if (haspackage==FALSE){
  install.packages("R2jags")
}
library("R2jags")


haspackage = require("nlme")
if (haspackage==FALSE){
  install.packages("nlme")
}
library("nlme")

# IQ Data Examples -----------------------------------------------------------------------------------------


# estimation with Least Squares:
model1LS = lm(y~1, data = dataIQ)
summary(model1LS)

# estimation with REML:
model1REML = gls(y~1, data = dataIQ, method= "REML")
summary(model1REML)
summary(model1REML)$sigma^2


# estimation with ML:
model1ML = gls(y~1, data = dataIQ, method= "ML")
summary(model1ML)
summary(model1ML)$sigma^2

# estimation with Bayesian

model01Bayes = function(){
  
  # likelihood
  for (i in 1:n){
    y[i] ~ dnorm(mu, tau)
  }
  
  #priors 
  mu ~ dnorm(100, 1/2.13^2)
  tau ~ dunif(1/400, 1/200)
  sigma2 = 1/tau
}

data = list(y = dataIQ$y, n = nrow(dataIQ))

jags.param = c("mu", "tau", "sigma2")

fit <- jags.parallel(data=data, 
                       parameters.to.save=jags.param,
                       n.iter=50000, n.chains=2,n.thin=2,n.burnin=40000,
                       model.file=model01Bayes)
chain = as.mcmc(fit)
plot(chain)
fit

model02Bayes = function(){
  
  # likelihood
  for (i in 1:n){
    y[i] ~ dnorm(mu, tau)
  }
  
  #priors 
  mu ~ dunif(-10000, 10000)
  tau ~ dunif(1/5000, 100000)
  sigma2 = 1/tau
}

data = list(y = dataIQ$y, n = nrow(dataIQ))

jags.param = c("mu", "tau", "sigma2")

fit2 <- jags.parallel(data=data, 
                     parameters.to.save=jags.param,
                     n.iter=50000, n.chains=2,n.thin=2,n.burnin=40000,
                     model.file=model02Bayes)
chain2 = as.mcmc(fit2)
plot(chain2)
fit2

model03Bayes = function(){
  
  # likelihood
  for (i in 1:n){
    y[i] ~ dnorm(mu, tau)
  }
  
  #priors 
  mu ~ dunif(-10000, 10000)
  tau ~ dgamma(.01, .01)
  sigma2 = 1/tau
}

data = list(y = dataIQ$y, n = nrow(dataIQ))

jags.param = c("mu", "tau", "sigma2")

fit3 <- jags.parallel(data=data, 
                      parameters.to.save=jags.param,
                      n.iter=50000, n.chains=2,n.thin=2,n.burnin=40000,
                      model.file=model03Bayes)
chain3 = as.mcmc(fit3)
plot(chain3)
fit3



model04Bayes = function(){
  
  # likelihood
  for (i in 1:n){
    y[i] ~ dnorm(mu, tau)
  }
  
  #priors 
  mu ~ dunif(102, 103)
  tau ~ dunif(1/242, 1/238)
  sigma2 = 1/tau
}

data = list(y = dataIQ$y, n = nrow(dataIQ))

jags.param = c("mu", "tau", "sigma2")

fit4 <- jags.parallel(data=data, 
                      parameters.to.save=jags.param,
                      n.iter=50000, n.chains=2,n.thin=2,n.burnin=40000,
                      model.file=model04Bayes)
chain4 = as.mcmc(fit4)
plot(chain4)
fit4
