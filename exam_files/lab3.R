set.seed(12345)
setwd('C:\\Users\\ekblo\\LiU\\tdde07')
library(coda)

#1
data = readRDS("Precipitation.rds")
log_data = log(data)
nDraws = 1000

#a

gibbsDraws = matrix(0, nDraws, 2)
mu_0 = 5
sigma_sqr_0 = 10
tau_sqr_0 = 5 #if tau big less sure about prior
n = length(data)
v_0 = 4

sigma_sqr = sigma_sqr_0
gibbs_matrix = matrix(0, nDraws, 2)

#Gibbs sampling
for (i in 1:nDraws) {
  #necessary parameters to draw mu, formulas from page 28
    w = (n / sigma_sqr) / ((n / sigma_sqr) + 1 / tau_sqr_0)
    mu_n = w * mean(log_data) + (1 - w) * mu_0
    tau_sqr_n = 1 / (n / sigma_sqr + 1 / tau_sqr_0)
    
    #draw mu (theta1)
    mu = rnorm(1, mean = mu_n, sd = tau_sqr_n)
    gibbs_matrix[i, 1] = mu
    
    #necessary parameters to draw sigma_sqr, formulas from page 44
    v_n = v_0 + n
    
    #draw sigma_sqr (theta2)
    sigma_sqr = v_n * ((v_0 * sigma_sqr_0 + sum((log_data - mu)^2)) / (n + v_0)) / rchisq(1, v_n)
    gibbs_matrix[i, 2] = sigma_sqr
    
}


#inefficiency factor

plot(gibbs_matrix[,1], main="mu", type="l")
plot(gibbs_matrix[,2], main="sigma_sqr", type="l")
effectiveSize(gibbs_matrix[,1])
effectiveSize(gibbs_matrix[,2])

auto_mu = acf(gibbs_matrix[,1])
auto_sigma_sqr = acf(gibbs_matrix[,2])
1 + 2 * sum(auto_mu$acf[-1])
1 + 2 * sum(auto_sigma_sqr$acf[-1])

#b
simulated_draws = matrix(0,nDraws-100, 1)
for (i in 100:nDraws) {
  simulated_draws[i] = rnorm(1, gibbs_matrix[i,1], gibbs_matrix[i,2])
}

plot(density(exp(simulated_draws)), col = "red", xlim = c(-5, 60))
lines(density(data))

#2
ebay_data = read.csv("eBayNumberOfBidderData.dat", header = TRUE, sep = "")

#a
beta_model = glm(formula=nBids ~ PowerSeller + VerifyID + Sealed + Minblem + MajBlem + LargNeg + LogBook + MinBidShare, data=ebay_data, family=poisson())
print(beta_model)
#b
library(mvtnorm)
#prior draws
X = as.matrix(ebay_data[,-1])
y = ebay_data$nBids
Sigma = 100 * solve(t(X)%*%X)
mu = matrix(0, 9, 1)
initVal = numeric(9)
  
  
LogPostPoisson <- function(betas,y,X,mu,Sigma){ #from lecture
  linPred <- (X%*%betas);
  logLik <- sum(y*linPred- exp(linPred))
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  
  return(logLik + logPrior)
  }
  
optimRes = optim(initVal,LogPostPoisson, gr=NULL, y, X, mu, Sigma, method=c("BFGS"), control=list(fnscale=-1), hessian=TRUE)

beta_tilde = optimRes$par
JInv = solve(-optimRes$hessian)

nDraws = 1000
betas_posterior = rmvnorm(nDraws, beta_tilde, JInv)
plot(density(betas_posterior[,9]))

#c

LogPostPoisson <- function(theta,y,X,mu,Sigma){ #from lecture
  linPred <- (X%*%t(theta));
  logLik <- sum(y*linPred- exp(linPred))
  logPrior <- dmvnorm(theta, mu, Sigma, log=TRUE);
  
  return(logLik + logPrior)
}

metropolisRandomWalk <- function(c, SIGMA, fcn, ...) { #... = y, X, mu
  n_param = as.numeric(dim(SIGMA)[1])
  theta = matrix(0,1,n_param);
  theta_matrix = matrix(0, n_param, nDraws)
  for (i in 1:nDraws) {
    theta_proposal = rmvnorm(1, theta, c*SIGMA) #step 1
    
    logPost_old = fcn(theta, ...) #step 2
    logPost_new = fcn(theta_proposal, ...) # pass them here
    alpha = min(1, exp(logPost_new - logPost_old))
    
    comparer = runif(1,0,1) #step 3
    if (alpha >= comparer) { 
      theta = theta_proposal
    }
    theta_matrix[,i] = theta
  }
  return (theta_matrix)
  
}
set.seed(12345)
nDraws = 1000
mu = matrix(0, 9, 1)
X = as.matrix(ebay_data[,-1])
y = ebay_data$nBids
Sigma = 100 * solve(t(X)%*%X)
betas = metropolisRandomWalk(1, JInv, LogPostPoisson, y, X, mu, Sigma)
plot(betas[9,], type="l") #plot 9th row (MinBidShare)

#d

Xd = matrix(c(1, 1, 0, 1, 0, 1, 0, 1.2, 0.8), 9, 1)
nBidders = c()
count = 0
for (i in 1:nDraws) {
  nBidders[i] = rpois(1, exp(Xd%*%t(betas[,i])))

}
length(nBidders[nBidders < 1]) / nDraws
max(nBidders)
plot(hist(nBidders, breaks=10))
library(ggplot2)
qplot(nBidders, geom="histogram")

#3
#a
mu = 13
sigma_sqr = 3
T = 300

ar_process <- function(phi) {
  X = c()
  X[1] = mu
  for (i in 2:T) {
    X[i] = mu + phi * (X[i-1] - mu) + rnorm(1, 0, sqrt(sigma_sqr))
  }
  return(X)
}

plot(ar_process(-1))
plot(ar_process(1))
plot(ar_process(0.1))
plot(ar_process(-0.0001))


#b
#i)
x = ar_process(0.2)
y = ar_process(0.95)

library(rstan)
#StanModel = stan_model('3d.stan')

StanModel = 'data {
  int<lower=0> T; //number of observations
  vector[T] y;
}

parameters {
  real<lower=-1, upper=1> phi;
  real mu;
  real<lower=0> sigma_sqr;
  
}
model {
  mu ~ normal(0,10); // Normal with mean 0, st.dev. 100
  sigma_sqr ~ scaled_inv_chi_square(1,2); // Scaled-inv-chi2 with nu 1,sigma 2
  for(i in 2:T){
    y[i] ~ normal(mu + phi*y[i-1], sqrt(sigma_sqr));
  }
}'

data_1 <- list(T=T, y=x)
warmup <- 1000
niter <- 2000
fit_x <- stan(model_code=StanModel, data=data_1, warmup=warmup,iter=niter,chains=4)
# Print the fitted model
print(fit_x,digits_summary=3)

data_2 <- list(T=T, y=y)
fit_y <- stan(model_code=StanModel, data=data_2, warmup=warmup,iter=niter,chains=4)
# Print the fitted model
print(fit_y,digits_summary=3)


#ii)
# Extract posterior samples
postDraws_x <- extract(fit_x)
postDraws_y <- extract(fit_y)


#plot joint posterior and convergence
plot(postDraws_x$mu, postDraws_x$phi)
plot(postDraws_x$mu, type='l')
plot(postDraws_x$phi, type='l')
plot(postDraws_x$sigma_sqr, type='l')

plot(postDraws_y$mu, postDraws_y$phi)
plot(postDraws_y$mu, type='l')
plot(postDraws_y$phi, type='l')
plot(postDraws_y$sigma_sqr,type='l')
