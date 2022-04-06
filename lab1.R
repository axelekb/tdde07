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
