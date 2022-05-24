#1 a)
data = read.csv("\\\\ad.liu.se\\home\\axeek668\\TDDE07\\exam_files\\TempLambohov.txt", header = TRUE, sep = "")
set.seed(12345)
library(mvtnorm)

#prior hyperparameters
mu0 = c(-10, 100, -100)
omega0 = matrix(0.8 * diag(3), 3, 3)
v0 = 3
sigma_sqr0 = 2

nDraws = 10
n = as.numeric(dim(data)[1])

#draws for prior sigma square ~ inverse chi-squared(v0, sigma_sqr0)
x = rchisq(nDraws, v0)
sigma_sqr = matrix(v0 * sigma_sqr0 / x)

#draws for prior betas ~ N(mu0, sigma_sqr*omega0^-1)
#beta_prior_draw = rmvnorm(nDraws, mu0, sigma_sqr0 * solve(omega0))

#for each prior draw of sigma square, draw prior betas ~ N(mu0, sigma_sqr * omega0^-1)
beta_prior_draw = apply(sigma_sqr, 1, function(x) rmvnorm(1, mu0, x * solve(omega0)))

#computing regression curve for each draw
prior_temp = matrix(nrow=nDraws, ncol=365)
for (i in 1:nDraws) {
  prior_temp[i,] = c(beta_prior_draw[1,i] + beta_prior_draw[2,i] * data$time + beta_prior_draw[3,i] * data$time^2)
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
joint_posterior_draw = function (nDraws, X, y){
  
  #setting up needed parameters
  beta_hat = solve(t(X) %*% X) %*% t(X) %*% y
  mun = (solve(t(X) %*% X + omega0)) %*% ((t(X) %*% X %*% beta_hat) + omega0 %*% mu0)
  omega_n = t(X) %*% X + omega0
  vn = v0 + n
  vn_sigma_sqr_n = v0 * sigma_sqr0 + (t(y) %*% y + t(mu0) %*% omega0 %*% mu0 - t(mun) %*% omega_n %*% mun)
  sigma_sqr_n = vn_sigma_sqr_n / vn
  
  #draw sigma_sqr ~ inverse chi-squared(vn, sigma_sqr_n)
  x = rchisq(nDraws, vn)
  sigma_sqr = as.numeric(vn * sigma_sqr_n) / x
  
  #draw betas ~ N(mun, sigma_sqr * omega_n^-1)
  beta = matrix(nrow = nDraws, ncol = 3)
  for (i in 1:nDraws) {
    beta[i,] = c(rmvnorm(1, mun, sigma_sqr[i] * solve(omega_n)))
  }
  
  draws = list(sigma_sqr, beta)
  names(draws) = c("sigma_sqr", "betas")
  return(draws)
}
  
X = matrix(nrow = 365, ncol = 3)
X[,1] = 1
X[,2] = data$time
X[,3] = data$time^2
draws = joint_posterior_draw(1000, X, data$temp)

# i
# inference of the posterior betas and sigma_sqr
plot(hist(draws$sigma_sqr))
plot(hist(draws$betas[,1]), xlab='beta0')
plot(hist(draws$betas[,2]), xlab='beta1')
plot(hist(draws$betas[,3]), xlab='beta2')

# ii
#the posterior regression function, each row is a day and each column a draw
f_time = X %*% t(draws$betas)
#taking median for each day (row) of the regression function
medians = apply(f_time, 1, median)
#95% equal tail intervals of the posterior regression function
bottom_95 = apply(f_time, 1, function(x) quantile(x,0.025))
top_95 = apply(f_time, 1, function(x) quantile(x, 0.975))

plot(data)
lines(data$time,medians, col="red")
lines(data$time, bottom_95, col="blue")
lines(data$time, top_95, col="blue")

#1c
#derivate of f_time w.r.t. time: B1 + 2B2*time = 0 => time = -B1/2B2 = x_tilde
#using this with the posterior draws of B1 and B2 to compute the posterior draws of x_tilde
x_tilde_draws = -draws$betas[,2] / (2*draws$betas[,3])
plot(density(365 * x_tilde_draws), xlab='x_tilde * 365', main='x_tilde posterior distribution')

#1d
# by using smoothness/shrinkage/regularization prior
# betas|sigma_sqr ~ N(0,sigma_sqr/lambda), i.e., mu_0=0 & omega_0=lambda*I
# --------------------------------------------------------------------------------------------------------------
#2a
library(mvtnorm)
set.seed(12345)
data = read.csv("\\\\ad.liu.se\\home\\axeek668\\TDDE07\\exam_files\\WomenAtWork.dat", header = TRUE, sep = "")

#setting parameters
tao = 5
nDraws = 1000

X = as.matrix(data[,-1])
y = data$Work

Npar = as.numeric(dim(X)[2])
initVal = matrix(0, Npar, 1)

#function for returning something proportional to the log-posterior for the betas that are inputted
LogPostLogistic <- function(betas,y,X,prior_mean,prior_sigma){ 
  linPred <- X %*% betas;
  logLik <- sum(linPred*y - log(1 + exp(linPred)));
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betas, prior_mean, prior_sigma, log=TRUE);
  return(logLik + logPrior)
}

#searching for beta_tilde and JInv by maximizing the posterior probability, with help of the function above
#using prior mean = 0 and prior sigma = tao^2*I for betas, according to instructions
optimRes = optim(initVal,LogPostLogistic, gr=NULL, y, X, as.matrix(rep(0,Npar)), tao^2*diag(Npar), method=c("BFGS"), control=list(fnscale=-1), hessian=TRUE)

JInv = solve(-optimRes$hessian)
beta_tilde = optimRes$par

#making the posterior draws
betas_posterior = rmvnorm(nDraws, beta_tilde, JInv)

#95%  equal tail posterior probability interval for NSmallChild
bottom_95 = quantile(betas_posterior[,6],0.025)
top_95 = quantile(betas_posterior[,6],0.975)

#maximum likelihood estimate for comparison
ml_estimate_model = glm(Work~0+., data = data, family = binomial)
ml_estimate_model$coefficients["NSmallChild"]

#b
set.seed(12345)
X = c(1, 20, 12, 8, 43, 0, 2) #woman described in lab

posterior_draws <- function(X, nDraws, beta_tilde, JInv) {
  #using the approximated normal posterior for beta done in a)
  betas_posterior = rmvnorm(nDraws, beta_tilde, JInv)
  #the logistic regression model, with the posterior betas and X plugged in
  posterior_predictive = exp(t(X) %*% t(betas_posterior)) / (1 + exp(t(X) %*% t(betas_posterior)))
  return(posterior_predictive)
}

woman_posterior_draws = posterior_draws(X, nDraws, beta_tilde, JInv)
plot(density(woman_posterior_draws))

#c
set.seed(12345)
women_working <- function(X, nDraws, beta_tilde, JInv) {
  betas_posterior = rmvnorm(nDraws, beta_tilde, JInv)
  posterior_predictive = exp(t(X) %*% t(betas_posterior)) / (1 + exp(t(X) %*% t(betas_posterior)))
  #posterior_mean = mean(posterior_predictive)  #using the posterior mean as estimated probability that the woman works?
  women_working = rbinom(nDraws, 11, posterior_predictive)
  return(women_working)
}
women_working = women_working(X, nDraws, beta_tilde, JInv)
plot(hist(women_working))
