set.seed(12345)
setwd("/home/gusli281/Documents/tdde07/lab3/")
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

