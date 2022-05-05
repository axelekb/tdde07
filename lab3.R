set.seed(12345)
setwd('C:/Users/Gustaf/OneDrive/Dokument/tdde07')
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
  theta = matrix(0,1,9);
  theta_matrix = matrix(0, 9, nDraws)
  for (i in 1:nDraws) {
    theta_proposal = rmvnorm(1, theta, c*SIGMA) #step 1
    
    logPost_old = fcn(theta, ...) #step 2
    logPost_new = fcn(theta_proposal, ...) # pass them here
    alpha = min(1, exp(logPost_old - logPost_new))
    
    comparer = runif(1,0,1) #step 3
    if (alpha < comparer) { 
      theta = theta_proposal
    }
    theta_matrix[,i] = theta
  }
  return (theta_matrix)
  
}
mu = matrix(0, 9, 1)
Sigma = 100 * solve(t(X)%*%X)
betas = metropolisRandomWalk(1, JInv, LogPostPoisson, y, X, mu, Sigma)
plot(betas[9,], type="l") #plot 9th row (MinBidShare)
