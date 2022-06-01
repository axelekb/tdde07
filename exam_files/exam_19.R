# 1a)
nDraws = 10000
posterior_4 = rbeta(nDraws, sqrt(4), 20)
posterior_16 = rbeta(nDraws, sqrt(16),20)
plot(density(posterior_4))
lines(density(posterior_16), col="red")

# 1b)
mean(posterior_4>=0.1)
mean(posterior_16>=0.1)

# 1c)
u = function(theta, c){
  return(100 + 20*log(theta) - c)
}
c_grid = seq(4,20,0.5)
utils = c()
for(i in 1:length(c_grid)){
  posterior = rbeta(nDraws, sqrt(c_grid[i]), 20)
  utils[i] = mean(sapply(posterior, u, c_grid[i]))
}
plot(c_grid, utils, type='l')
points(c_grid[which.max(utils)],utils[which.max(utils)], col="red")
c_grid[which.max(utils)]

# 2a)
log_posterior = function(theta, n, x) {
  log_likelihood = sum(log(dbinom(x, n, theta)))
  #log_likelihood = sum(x)*log(theta) + sum(n-x)*log(1-theta)
  log_prior = 2*log(1-theta)

  return(log_likelihood + log_prior)
}
theta_grid = seq(0,1,0.001)
posterior = c()
for(i in 1:length(theta_grid)){
  posterior[i] = log_posterior(theta_grid[i], 50, ebay)
}

plot(theta_grid, exp(posterior)/(0.001*sum(exp(posterior))), type='l')

# 2b)
set.seed(100)
K = 2
nIter = 500
gibbs_sample = GibbsMixPois(ebay, 2, 1, 1, 1, xGrid, 500)
plot(gibbs_sample$thetaSample[,1], type='l')
plot(gibbs_sample$thetaSample[,2], type='l')
plot(cumsum(gibbs_sample$thetaSample[,1])/seq(1:length(gibbs_sample$thetaSample[,1])), type='l')
plot(cumsum(gibbs_sample$thetaSample[,2])/seq(1:length(gibbs_sample$thetaSample[,2])), type='l')

