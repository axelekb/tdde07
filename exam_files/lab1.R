# 1.a)
s = 13
n = 50
f = n - s
a0 = 5
b0 = 5

#the different amount of draws
draws = matrix(seq(from=10, to= 10000, by=10))

#function for making the draws and calculate mean and sd of theta
draw = function(nDraws, s, n, f, a0, b0){
  posterior = rbeta(nDraws, a0 + s, b0 + f)
  mean = mean(posterior)
  sd = sd(posterior)
  return(c(mean, sd))
}

#for each amount of draws, calculate means and sd's
means_sd = apply(draws, 1, draw, s=s, n=n, f=f, a0=a0, b0=b0)

#plot results
plot(means_sd[1,])
plot(means_sd[2,])


# 1.b)
nDraws = 10000

#make the draws
x = rbeta(nDraws, a0 + s, b0 + f)

#check fraction of the draws smaller than 0.3, becomes the posterior probability of theta smaller than 0.3
sum(x < 0.3) / length(x)

#exact value from the beta cdf
pbeta(0.3, a0 + s, b0 + f)

# 1.c)
nDraws = 10000
x = rbeta(nDraws, a0 + s, b0 + f)
x = log(x/(1-x))
plot(density(x))

# --------------------------------------------------------------------------------------------------------------

# 2.a)
nDraws = 10000
mu = 3.5
sample = c(33, 24, 48, 32, 55, 74, 23, 76, 17)
n = as.numeric(length(sample))

#drawing sigma_sq (posterior) from inv-chi(n,tau_sq)
tau_sq = sum((log(sample) - mu)^2)/n
x = rchisq(nDraws, n)
sigma_sq = n * tau_sq / x

plot(density(sigma_sq))

# 2.b)
# calculating gini-coefficient
G = 2 * pnorm(sqrt(sigma_sq / 2)) - 1
plot(density(G))

# 2.c)
#computing an empirical cdf for the posterior draws, and with the help of this printing the 2.5% and 97.5% quantiles for the posterior density(G)
cdf = ecdf(G)
quantile(cdf,0.025)
quantile(cdf,0.975)

# 2.d)
density = cbind(density(G)$x, density(G)$y)
density = density[order(density[,2],decreasing=TRUE),]
density = cbind(density, cumsum(density[,2])/sum(density[,2]))

#density now contains in decreasing order of density: parameter value, density, cumulative density
#now selecting the parameter values (column 1) where the proportional cumulative sum of densities (column 3) is less than 95%
highest_densites = density[1:sum(density[,3]<0.95),1]

#the min/max parameter value of the 95% highest density values is the left/right HPDI bound
min(highest_densites)
max(highest_densites)

# 3.a)
# p(??|y, ?) ??? p(y|?, ??) * p(??)
x = c(1.83, 2.02, 2.33, -2.79, 2.07, 2.02, -2.44, 2.14, 2.54, 2.23)
lambda = 1
kappas = seq(0.01, 10, 0.01)
mu = 2.51

likelihood = c()
prior = c()
for (kappa in seq(0.01, 10, 0.01)){
  likelihood = append(likelihood, prod(exp(kappa * cos(x - mu))/(2 * pi * besselI(kappa, 0))))
  prior = append(prior, lambda * exp(-lambda * kappa))
}
posterior = likelihood * prior
plot(kappas,posterior/(sum(posterior)*0.01), ylab="Posterior probability", xlab="Kappa", type="l")

# 3.b)
kappas[which.max(posterior)]
