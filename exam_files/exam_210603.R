#exam 21-06-03
# 1a)
n = 100
s = 38
f = 100 - s
alfa = 16
beta = 24

nDraws = 1000
posterior_draws = rbeta(nDraws, alfa + s, beta + f)

1 - pbeta(0.4, alfa + s, beta + f)

plot(density((1 - posterior_draws)/ posterior_draws))

# 1b)
ratio = (1 - posterior_draws) / posterior_draws
quantile(ratio, 0.025)
quantile(ratio, 0.975)

# 1c)
beta(alfa + s, beta + f) / beta(alfa, beta)

# 1d)
set.seed(12345)
alpha_1 = 20
alpha_2 = 20
alpha_3 = 20

y1 = 38
y2 = 27
y3 = 35

nDraws = 1000
posterior = matrix(nrow=nDraws, ncol=3)
posterior[,1] = rgamma(nDraws, alpha_1 + y1, 1)
posterior[,2] = rgamma(nDraws, alpha_2 + y2, 1)
posterior[,3] = rgamma(nDraws, alpha_3 + y3, 1)
for (i in 1:nDraws){
  sum = sum(posterior[i,])
  posterior[i,1] = posterior[i,1] / sum
  posterior[i,2] = posterior[i,2] / sum
  posterior[i,3] = posterior[i,3] / sum
}
mean(posterior[,1] > posterior[,3])

# -------------------------------------------------------------------------------------
# 3a)
nDraws = 10000
mu0 = rep(0, 7)
omega0 = 1/25 * diag(7)
v0 = 1
sigma_sqr0 = 2^2

posterior_draws = BayesLinReg(y, X, mu0, omega0, v0, sigma_sqr0, nDraws)

apply(posterior_draws$betaSample, 2, mean)
apply(posterior_draws$betaSample, 2, function(x) quantile(x, c(0.025, 0.975)))
quantile(posterior_draws$betaSample[,2], c(0.025, 0.975))
# the regression coefficient for verbal IQ is with 95% probability between 0.5 and 0.9, meaning it is of importance for the test result

# 3b)
median(sqrt(posterior_draws$sigma2Sample))

# 3c)
plot(density(posterior_draws$betaSample[,6]-posterior_draws$betaSample[,7]))
# majority of probability mass above 0, meaning it is probable that the effect on y from x1 is different for
# students in high school B compared to C.

# 3d)
x1_grid = seq(min(X[,2]), max(X[,2]), by=0.01)
betas = BayesLinReg(y, X, mu0, omega0, v0, sigma_sqr0, nDraws)$betaSample
draws = matrix(NA, nDraws, length(x1_grid))
for (x1 in 1:length(x1_grid)){
  x = c(1, x1_grid[x1], 0.5, 0, 0, x1_grid[x1]*0, x1_grid[x1]*0)
  draws[,x1] = betas %*% x
}
intervals = apply(draws, 2, function(x) quantile(x, c(0.05, 0.95)))
plot(x1_grid, intervals[2,], type='l', ylim=c(min(intervals[1,]), max(intervals[2,])))
lines(x1_grid, intervals[1,])

# 3e)
x = c(1, 0.4, 1, 1, 0, 0.4*1, 0.4*0)
predictive_draws = betas %*% x + rnorm(10, 0, sqrt(posterior_draws$sigma2Sample))
plot(density(predictive_draws))
