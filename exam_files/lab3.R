set.seed(12345)
setwd('C:\\Users\\ekblo\\LiU\\tdde07')

#1
data = readRDS("\\\\ad.liu.se\\home\\axeek668\\TDDE07\\exam_files\\Precipitation.rds")
log_data = log(data)
nDraws = 1000

#a
mu_0 = 1
sigma_sqr_0 = 10
tau_sqr_0 = 5 #if tau big less sure about prior
n = length(data)
v_0 = 4

# initialising sigma_sqr
sigma_sqr = 1
gibbs_matrix = matrix(0, nDraws, 2)

#Gibbs sampling
for (i in 1:nDraws) {
    #necessary parameters to draw mu (as in case with normal data, normal prior and known variance since we assume theta2 (sigma_sqr) to be known)
    w = (n / sigma_sqr) / (n / sigma_sqr + 1 / tau_sqr_0)
    mu_n = w * mean(log_data) + (1 - w) * mu_0
    tau_sqr_n = 1 / (n / sigma_sqr + 1 / tau_sqr_0)
    
    #draw mu (theta1)
    mu = rnorm(1, mean = mu_n, sd = tau_sqr_n)
    gibbs_matrix[i, 1] = mu
    
    #necessary parameters to draw sigma_sqr
    v_n = v_0 + n
    
    #draw sigma_sqr (theta2)
    sigma_sqr = v_n * ((v_0 * sigma_sqr_0 + sum((log_data - mu)^2)) / (n + v_0)) / rchisq(1, v_n)
    gibbs_matrix[i, 2] = sigma_sqr
    
}

plot(gibbs_matrix[,1], main="mu", type="l")
plot(gibbs_matrix[,2], main="sigma_sqr", type="l")


#effective sample size and inefficiency factor 
library(coda)
effectiveSize(gibbs_matrix[,1])
effectiveSize(gibbs_matrix[,2])

auto_mu = acf(gibbs_matrix[,1])
auto_sigma_sqr = acf(gibbs_matrix[,2])
1 + 2 * sum(auto_mu$acf[-1])
1 + 2 * sum(auto_sigma_sqr$acf[-1])

#b
plot(density(data))

#predictive draws, one draw for each draw of mu and sigma_sqr from the gibbs sampling
cut = 100 # cut away the first 100 simulations from gibbs due to late convergence (warmup)
simulated_draws = matrix(0,nDraws-cut, 1)
for (i in 1:nDraws-cut) {
  simulated_draws[i] = rnorm(1, gibbs_matrix[i+cut,1], gibbs_matrix[i+cut,2])
}
lines(density(exp(simulated_draws)), col = "red")

# ------------------------------------------------------------------------------------------------------ 
#2
ebay_data = read.csv("\\\\ad.liu.se\\home\\axeek668\\TDDE07\\exam_files\\eBayNumberOfBidderData.dat", header = TRUE, sep = "")

#a
linear_model = glm(formula=nBids ~ PowerSeller + VerifyID + Sealed + Minblem + MajBlem + LargNeg + LogBook + MinBidShare, data=ebay_data, family=poisson())
summary(linear_model)
# coefficients with p-value below 0.05 are considered significant, there are 6 in this case
#b
library(mvtnorm)
#prior draws

  
LogPostPoisson <- function(betas,y,X,prior_mean,prior_sd){
  linPred <- (X%*%betas)
  logLik <- sum(y*linPred- exp(linPred)) # log-likelihood for the poisson regression
  logPrior <- dmvnorm(betas, prior_mean, prior_sd, log=TRUE) # normal prior
  return(logLik + logPrior)
}

X = as.matrix(ebay_data[,-1])
y = ebay_data$nBids

#prior parameters
Sigma = 100 * solve(t(X)%*%X)
mu = matrix(0, 9, 1)

initVal = numeric(9)
#maximising the posterior prob for the betas, to find the posterior mode and hessian
optimRes = optim(initVal,LogPostPoisson, gr=NULL, y, X, mu, Sigma, method=c("BFGS"), control=list(fnscale=-1), hessian=TRUE)

beta_tilde = optimRes$par
JInv = solve(-optimRes$hessian)

nDraws = 1000
#making the posterior beta draws from the approx. multivariate normal, using mode and JInv calculated above
betas_posterior = rmvnorm(nDraws, beta_tilde, JInv)
plot(density(betas_posterior[,9]))

#c
LogPostPoisson <- function(theta,y,X,prior_mean,prior_sd){
  linPred <- (X%*%t(theta));
  logLik <- sum(y*linPred- exp(linPred))
  logPrior <- dmvnorm(theta, prior_mean, prior_sd, log=TRUE);
  return(logLik + logPrior)
}

metropolisRandomWalk <- function(c, proposal_sigma, fcn, ...) { #... = y, X, prior mean, prior sigma
  n_param = as.numeric(dim(proposal_sigma)[1])
  theta = matrix(0,1,n_param);
  theta_matrix = matrix(0, n_param, nDraws)
  for (i in 1:nDraws) {
    theta_proposal = rmvnorm(1, theta, c*proposal_sigma) #sample proposal
    logPost_old = fcn(theta, ...) #posterior probability for old theta
    logPost_new = fcn(theta_proposal, ...) #posterior probability for proposal theta
    alpha = min(1, exp(logPost_new - logPost_old)) #acceptance probability
    
    comparer = runif(1,0,1) #with probability alpha, set theta_i = theta_proposal, theta_i-1 otherwise (i.e. no update)
    if (alpha >= comparer) { 
      theta = theta_proposal
    }
    theta_matrix[,i] = theta
  }
  return (theta_matrix)
  
}
set.seed(12345)
nDraws = 1000
X = as.matrix(ebay_data[,-1])
y = ebay_data$nBids

prior_mean = matrix(0, 9, 1)
prior_sd = 100 * solve(t(X)%*%X)

c = 1
betas = metropolisRandomWalk(c, JInv, LogPostPoisson, y, X, prior_mean, prior_sd)

par(mfrow = c(3,3))
plot(betas[1,], type="l")
plot(betas[2,], type="l")
plot(betas[3,], type="l")
plot(betas[4,], type="l")
plot(betas[5,], type="l")
plot(betas[6,], type="l")
plot(betas[7,], type="l")
plot(betas[8,], type="l")
plot(betas[9,], type="l") 
#after roughly 400 iterations, all parameters have converged

#d
X = matrix(c(1, 1, 0, 1, 0, 1, 0, 1.2, 0.8), 1, 9)
nBidders = c()
cut = 400
for (i in 1:nDraws-cut) {
  nBidders[i] = rpois(1, exp(t(X) %*% betas[,i + cut])) #for each draw of betas, draw predictive y_new given those betas
}
length(nBidders[nBidders == 0]) / nDraws #using the predictive draws to calculate p(y_new=0|y)

dev.off()
hist(nBidders, right="F") # histogram showing wrong values on y-axis??
library(ggplot2)
qplot(nBidders, geom="histogram") # this is correct

# -----------------------------------------------------------------------------------------------------------
#3
#a
mu = 13
sigma_sqr = 3
T = 300

ar_1_process <- function(phi) {
  X = c()
  X[1] = mu
  for (i in 2:T) {
    X[i] = mu + phi * (X[i-1] - mu) + rnorm(1, 0, sqrt(sigma_sqr))
  }
  return(X)
}

plot(ar_1_process(-1))
plot(ar_1_process(1))
plot(ar_1_process(0.1))
plot(ar_1_process(-0.0001))


#b
#i)
x = ar_1_process(0.2)
y = ar_1_process(0.95)

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
    y[i] ~ normal(mu + phi*(y[i-1]-mu), sqrt(sigma_sqr));
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

#in general, smaller variance in estimated parameters when phi=0.2 and following this, a parameter estimation that works better in general 

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
