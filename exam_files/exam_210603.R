#exam 21-06-03
# 1a)
n = 100
s = 38
f = 100 - s
alfa = 16
beta = 24

nDraws = 1000
posterior_draws = rbeta(nDraws, alfa + s, beta + f)

1 - pbeta(0.4, alfa + s, beta + f)

plot(density((1 - posterior_draws)/ posterior_draws))

# 1b)
ratio = (1 - posterior_draws) / posterior_draws
quantile(ratio, 0.025)
quantile(ratio, 0.975)

# 1c)
beta(alfa + s, beta + f) / beta(alfa, beta)

# 1d)
set.seed(12345)
alpha_1 = 20
alpha_2 = 20
alpha_3 = 20

y1 = 38
y2 = 27
y3 = 35

nDraws = 1000
posterior = matrix(nrow=nDraws, ncol=3)
posterior[,1] = rgamma(nDraws, alpha_1 + y1, 1)
posterior[,2] = rgamma(nDraws, alpha_2 + y2, 1)
posterior[,3] = rgamma(nDraws, alpha_3 + y3, 1)
for (i in 1:nDraws){
  sum = sum(posterior[i,])
  posterior[i,1] = posterior[i,1] / sum
  posterior[i,2] = posterior[i,2] / sum
  posterior[i,3] = posterior[i,3] / sum
}
mean(posterior[,1] > posterior[,3])

# -------------------------------------------------------------------------------------
# 2
