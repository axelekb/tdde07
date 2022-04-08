# 1.a)
s = 13
n = 50
f = n - s
a0 = 5
b0 = 5

means = c()
deviations = c()
for (n in seq(from=10, to= 10000, by=10)){
  set.seed(12345)
  posterior = rbeta(n, a0 + s, b0 + f)
  means <- append(means, mean(posterior))
  deviations <- append(deviations, sd(posterior))
}
plot(means)
plot(deviations)

# 1.b)
nDraws = 10000
x = rbeta(nDraws, a0 + s, b0 + f)
t = 0
for (n in 1:length(x)) {
  if (x[n] < 0.3) {
    t = t + 1
  }
}
t/length(x)
pbeta(0.3, a0 + s, b0 + f)

# 1.c)
set.seed(12345)
nDraws = 10000
x = rbeta(nDraws, a0 + s, b0 + f)
x = log(x/(1-x))
plot(density(x))

# -------------------------------------------------------

# 2.a)
set.seed(12345)
nDraws = 10000
mean = 3.5
sample = c(33, 24, 48, 32, 55, 74, 23, 76, 17)
n = as.numeric(length(sample))
tau_sqr = sum((log(sample) - mean)^2)/n

x = rchisq(nDraws, n)
sigma_sq = n * tau_sqr / x
plot(density(sigma_sq))

# 2.b)
G = 2*pnorm(sqrt(sigma_sq/2)) - 1


# 2.c)
cdf = ecdf(G)
quantile(cdf,0.025)
quantile(cdf,0.975)

# 2.d)
plot(density(G))
density = density(G)
densities = sort(density(G)[["y"]],decreasing=TRUE)

density_sum = sum(densities)
cumulative_sum = cumsum(densities)
# searching for the index of the values such that the cumulative sum is less than 95% of total sum
h_index = order(density$y,decreasing=TRUE)[1:length(cumulative_sum[cumulative_sum < 0.95 * density_sum])]

min(density$x[h_index])
max(density$x[h_index])


# 3.a)
# p(??|y, µ) ??? p(y|µ, ??) * p(??)
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
plot(posterior)
plot(kappas,posterior/(sum(posterior)*0.01))

# 3.b)
kappas[which.max(posterior)]
