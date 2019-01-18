irtItemProb = function(a, b, c=0, theta){
  prob = c + (1-c) * exp(a*(theta-b))/(1+exp(a*(theta-b)))
  return(prob)
}

trueTheta = 0
nItems = 15
bRange = c(-2,2)
aRange = c(1,2)
bSE = 1
aSE = 1
nSamples = 100000

# draw mean values of a, b
a = runif(n = nItems, min = aRange[1], max = aRange[2])
b = runif(n = nItems, min = bRange[1], max = bRange[2])

# draw items
itemResponses = rbinom(n = nItems, size = 1, prob = irtItemProb(a = a, b = b, theta = 1))
thetaChain = list(rep(NA, nSamples), rep(NA, nSamples))

# initialize theta values
curTheta = trueTheta
curThetaRand = trueTheta
for (iteration in 1:nSamples){
  
  # draw item parameters (if random)
  iterA = rnorm(n = nItems, mean = a, sd = aSE)
  iterB = rnorm(n = nItems, mean = b, sd = bSE)
  
  # calculate current likelihood of the data | theta
  curLogLike = sum(dbinom(x = itemResponses, size = 1, prob = irtItemProb(a = a, b = b, theta = curTheta), log = TRUE))
  curLogLikeRand = sum(dbinom(x = itemResponses, size = 1, prob = irtItemProb(a = iterA, b = iterB, theta = curThetaRand), log = TRUE))
  
  # draw new theta value
  propTheta = rnorm(n = 1, mean = curTheta, sd = 1)
  propThetaRand = rnorm(n = 1, mean = curThetaRand, sd = 1)
  
  # calculate proposed likelihood of the data | theta
  propLogLike = sum(dbinom(x = itemResponses, size = 1, prob = irtItemProb(a = a, b = b, theta = propTheta), log = TRUE))
  propLogLikeRand = sum(dbinom(x = itemResponses, size = 1, prob = irtItemProb(a = iterA, b = iterB, theta = propThetaRand), log = TRUE))
  
  # do MH:
  if (log(runif(n = 1)) < (propLogLike-curLogLike)){
    # accept
    curTheta = propTheta
  } 

  # do MH:
  if (log(runif(n = 1)) < (propLogLikeRand-curLogLikeRand)){
    # accept
    curThetaRand = propThetaRand
  }
  
  thetaChain[[1]][iteration] = curTheta
  thetaChain[[2]][iteration] = curThetaRand
}

par(mfrow = c(1,2))

plot(thetaChain[[1]], type="l", ylab = expression(theta), xlab = "Iteration Number")
lines(thetaChain[[2]], type="l", col = 2)
plot(density(thetaChain[[1]]), col = 1, main="")
lines(density(thetaChain[[2]]), col = 2)

par(mfrow = c(3,2))
plot(thetaChain[[1]], type="l", ylab = expression(theta), xlab = "Iteration Number")
plot(thetaChain[[2]], type="l", ylab = expression(theta), xlab = "Iteration Number", col =2)
plot(density(thetaChain[[1]]), col = 1, main="")
plot(density(thetaChain[[2]]), col = 2, main="")
plot(thetaChain[[1]], type="l", ylab = expression(theta), xlab = "Iteration Number")
lines(thetaChain[[2]], type="l", col = 2)
plot(density(thetaChain[[1]]), col = 1, main="")
lines(density(thetaChain[[2]]), col = 2)
