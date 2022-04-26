#1 a)
set.seed(12345)
library(mvtnorm)
u0 = c(-10, 100, -100)
omega0 = matrix(0.02 * diag(3), 3, 3) #0.02 default
v0 = 3
sigma_sqr0 = 0.1
nDraws = 1000
data = read.csv("/home/gusli281/Documents/tdde07/TempLambohov.txt", header = TRUE, sep = "")
n = as.numeric(dim(data)[1])

x = rchisq(nDraws, v0)
#s_sqr = sum((x - mean(x))^2)/(v0)
sigma_sqr = v0*sigma_sqr0 / x
plot(density(sigma_sqr), xlim=c(0,20)) #prior reasonable, majority between 0 and 5, ask later
beta_prior_draw = rmvnorm(nDraws, u0, sigma_sqr0*solve(omega0))

prior_temp = matrix(nrow=nDraws, ncol=365)
for (i in 1:nDraws) {
  prior_temp[i,] = c(beta_prior_draw[i,1] + beta_prior_draw[i,2] * data$time + beta_prior_draw[i,3] * data$time^2)
}
plot(data$time, prior_temp[1,], type='l', ylim=range(prior_temp))
lines(data$time,prior_temp[2,],col='red')
lines(data$time,prior_temp[3,],col='blue')
lines(data$time,prior_temp[4,],col='green')
lines(data$time,prior_temp[5,],col='purple')
lines(data$time,prior_temp[6,],col='yellow')
lines(data$time,prior_temp[7,],col='orange')
lines(data$time,prior_temp[8,],col='pink')
lines(data$time,prior_temp[9,],col='light blue')
lines(data$time,prior_temp[10,],col='grey')

#1b
sigma_sqr_posterior_draw = function (nDraws, X, y){

  y = data$temp
  beta_hat = solve(t(X) %*% X) %*% t(X) %*% y
  un = (solve(t(X) %*% X + omega0)) %*% ((t(X) %*% X %*% beta_hat) + omega0 %*% u0)
  omega_n = t(X) %*% X + omega0
  vn = v0 + n
  vn_sigma_sqr_n = v0 * sigma_sqr0 + (t(y) %*% y + t(u0) %*% omega0 %*% u0 - t(un) %*% omega_n %*% un)
  sigma_sqr_n = vn_sigma_sqr_n / vn
  
  x = rchisq(nDraws, vn)
  sigma_sqr = as.numeric(vn_sigma_sqr_n) / x
  
  beta = matrix(nrow = nDraws, ncol = 3)
  for (i in 1:nDraws) {
    beta[i,] = c(rmvnorm(1, un, sigma_sqr[i] * solve(omega_n)))
  }
  draws = list(sigma_sqr, beta)
  names(draws) = c("sigma_sqr", "betas")
  return(draws)
}
  
X = matrix(nrow = 365, ncol = 3)
X[,1] = 1
X[,2] = data$time
X[,3] = data$time^2
draws = sigma_sqr_posterior_draw(1000, X, data$temp)

# i
plot(hist(draws$sigma_sqr))
plot(hist(draws$betas[,1]), xlab='beta0')
plot(hist(draws$betas[,2]), xlab='beta1')
plot(hist(draws$betas[,3]), xlab='beta2')

# ii
f_time = X %*% t(draws$betas)
medians = apply(f_time, 1, median)
bottom_95 = apply(f_time, 1, function(x) quantile(x,0.025))
top_95 = apply(f_time, 1, function(x) quantile(x, 0.975))

plot(data)
lines(data$time,medians, col="red")
lines(data$time, bottom_95, col="blue")
lines(data$time, top_95, col="blue")

#1c
# derivate of f_time w.r.t. time: B1 + 2B2*time = 0 => time = -B1/2B2
plot(density(-365*draws$betas[,2]/(2*draws$betas[,3])), xlab='x_tilde * 365', main='x_tilde posterior distribution')

#2a
library(mvtnorm)
set.seed(12345)
data = read.csv("WomenAtWork.dat", header = TRUE, sep = "")
t = 5
Npar = 7 #matrix size
Sigma = t^2*diag(Npar)
X = as.matrix(data[,-1]) #all rows not first col
nDraws = 1000
mu <- as.matrix(rep(0,Npar))
initVal = matrix(0, Npar, 1)
y = data$Work

LogPostLogistic <- function(betas,y,X,mu,Sigma){ #from lecture
  linPred <- X%*%betas;
  logLik <- sum( linPred*y - log(1 + exp(linPred)) );
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  
  return(logLik + logPrior)
}

betas = rmvnorm(nDraws, mu, t*diag(Npar))
betaPrior = dmvnorm(X, mu, prior_sigma)
#par = posterior mode
optimRes = optim(initVal,LogPostLogistic, gr=NULL, y, X,mu, Sigma, method=c("BFGS"), control=list(fnscale=-1), hessian=TRUE)
#names(optimRes$par) = colnames(X)
JInv = solve(-optimRes$hessian)
beta_tilde = optimRes$par
betas_posterior = rmvnorm(nDraws, beta_tilde, JInv)

#95% for NSmallChild
bottom_95_child = quantile(betas_posterior[,6],0.025)
top_95_child = quantile(betas_posterior[,6],0.975)

glm.model = glm(Work~0+., data = data, family = binomial)

#b
set.seed(12345)
X_43 = c(1, 20, 12, 8, 43, 1, 1) #woman described in lab
probs_woman <- function(X, nDraws, beta_tilde, JInv) {
  betas_posterior = rmvnorm(nDraws, beta_tilde, JInv)
  probs = exp(t(X_43) %*% t(betas_posterior)) / (1 + exp(t(X_43) %*% t(betas_posterior)))
  return(probs)
}

probs = probs_woman(X_43, nDraws, beta_tilde, JInv)
plot(density(probs))

#c
set.seed(12345)
women_working <- function(X, nDraws, beta_tilde, JInv) {
  betas_posterior = rmvnorm(nDraws, beta_tilde, JInv)
  probs = exp(t(X_43) %*% t(betas_posterior)) / (1 + exp(t(X_43) %*% t(betas_posterior)))
  probs_avg = mean(probs)
  women_working = rbinom(nDraws, 11, probs_avg)
  return(women_working)
}
women_working = women_working(X, nDraws, beta_tilde, JInv)
plot(hist(women_working))
