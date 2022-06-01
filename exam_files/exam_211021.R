#  1d)
log_post = function(theta, n, xsum){
  return(log(theta^(2+xsum)*exp(-theta*(n+1/2))))
}
theta_grid = seq(2.5, 7.5, 0.01)

log_posterior_draws = log_post(theta_grid, 15, 75)
plot(theta_grid, exp(log_posterior_draws)/(0.01*sum(exp(log_posterior_draws))))

# 1e)
theta_init = 0
optimres = optim(4, log_post, gr=NULL, 15, 75, method=c("L-BFGS-B"), lower=3, control=list(fnscale=-1), hessian=TRUE)

theta_tilde = optimres$par
Jinv = solve(-optimres$hessian)
nDraws = 10000
posterior_approx = dnorm(theta_grid, theta_tilde, sqrt(Jinv))
lines(theta_grid, posterior_approx, col="red")

# 1f)
xsum = 75
n = 15
samples = matrix(NA, nDraws, n)
for (i in 1:nDraws){
  posterior_thetas = rgamma(n, 3 + xsum, n + 1/2)
  samples[i,] = poisson_draws = rpois(n, posterior_thetas)
}
T_xrep = apply(samples, 1, max)
mean(T_xrep >= 14)                             

# 2a)
nIter = 10000
posterior_draw = BayesLogitReg(y, X, rep(0,3), 16*diag(3), nIter)
betas = posterior_draw$betaSample
quantile(betas[,2], c(0.05, 0.95))
# 90% posterior probability that beta1 is on this interval

# 2b)
mean(betas[,3]>0)
# x2 is probable to have a positive effect on p since the regression coefficient to x2 is probably larger than 0

# 2c)
mean(betas[,2]>0 & betas[,3]>0)

# 2d)
xj = c(1, 0, 0)
patient_posterior = exp(betas %*% xj) / (1 + exp(betas %*% xj))
plot(density(patient_posterior))
mean(patient_posterior>0.5)

# 2e)
x1_grid = seq(min(X[,2]), max(X[,2]), 0.01)
intervals = matrix(NA, length(x1_grid), 2)
for(i in 1:length(x1_grid)){
  xk = c(1, x1_grid[i], 1)
  p_draws = exp(betas %*% xk) / (1 + exp(betas %*% xk))
  intervals[i,] = quantile(p_draws, c(0.025, 0.975))
}

plot(x1_grid, intervals[,2], ylim=c(min(intervals[,1]), max(intervals[,2])), type='l')
lines(x1_grid, intervals[,1])


