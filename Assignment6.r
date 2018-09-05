#Q1. RNormal

p = 0.05
m = 500
z = (rbinom(20000, size = m, prob = p) - m * p)/sqrt(m*p*(1-p))
qqnorm(z, ylim = c(-4, 4), main = paste('QQ-plot, m = ', m, 'p = ', p))
qqline(z)

qqnorm(z, main = paste('QQ-plot, m = ', m, 'p = ', p))
qqline(z)

par(mfrow = c(1,2))


for(m in 90:100){
  z = (rbinom(20000, size = m, prob = p)/sqrt(m*p*(1-p)))
  qqnorm(z, ylim = c(-4, 4), main = paste('QQ-plot, m = ', m, 'p = ', p))
  qqline(z)
}

#Q2. RPoisson
m = 119
z = (rpois(20000, lambda=m) - m)/sqrt(m)
qqnorm(z, ylim = c(-4, 4), main = 'QQ-plot')
qqline(z)
mtext(bquote(lambda == .(m)), 3)

for (m in seq(1, 120, 2))
{
  z = (rpois(20000, lambda=m) - m)/sqrt(m)
  qqnorm(z, ylim = c(-4, 4), main = 'QQ-plot')
  qqline(z)
  mtext(bquote(lambda == .(m)), 3)
}

#Q3. RCauchy
n = 2000
par(mfrow = c(2,2))
for (i in 1:4)
{	x = rcauchy(n)
x = cumsum(x)/1:n
plot(1:n, x, type = "l")
abline(h = 0, col = "red")
}
par(mfrow = c(1,1))
cauchy.mean = c()
for(i in 1:100){
  cauchy.mean = c(cauchy.mean, mean(rcauchy(n)))
}
hist(cauchy.mean, xlim = c(-5, 5), breaks = 100)



par(mfrow = c(1,1))
#Q5. Bootstrap & Jackknife
library(bootstrap)
data = c(57, 60, 52, 49, 56, 46, 51, 63, 49, 57)
theta1 = function(x) {mean(x)}
results.bootstrap = bootstrap(data, 1000, theta = theta1)
lowlimit = quantile(results.bootstrap$thetastar, 0.05)
upperlimit = quantile(results.bootstrap$thetastar, 0.95)
lowlimit
upperlimit


data = c(57, 60, 52, 49, 56, 46, 51, 63, 49, 57)
theta1 = function(x) {mean(x)}
results.jack = jackknife(data, theta1)
results.jack
results.jack$jack.se
lowlimit = mean(results.jack$jack.values) - 2 * results.jack$jack.se
upperlimit = mean(results.jack$jack.values) + 2 * results.jack$jack.se
lowlimit
upperlimit


# g(x) = C e^{-x^1.5}, x > 0, C: unknown
# rejection sampling: k=0.5/C
kg = function(x) 0.5 * exp(-(x^1.5))
X = rexp(100000)
U = runif(100000)
X = X[ U*dexp(X) < kg(X) ]
hist(X, prob=T, breaks="Scott")


## law of large numbers
n = 1000
x1 = runif(n)
x1 = cumsum(x1)/1:n
x2 = rnorm(n)
x2 = cumsum(x2)/1:n
x3 = rexp(n)
x3 = cumsum(x3)/1:n
x4 = rpois(n, 1)
x4 = cumsum(x4)/1:n

par(mfrow = c(2,2))
plot(1:n, x1, type = "l", main = "U(0,1)")
abline(h = 0.5, col = "red")
plot(1:n, x2, type = "l", main = "N(0,1)")
abline(h = 0, col = "red")
plot(1:n, x3, type = "l", main = "exp(1)")
abline(h = 1, col = "red")
plot(1:n, x4, type = "l", main = "Poisson(1)")
abline(h = 1, col = "red")

x = 1/sqrt(runif(n))
x = cumsum(x)/1:n
win.graph()
plot(1:n, x, type = "l", main = "Infinite variance")
abline(h = 2, col = "red")

par(mfrow = c(2,2))
for (i in 1:4)
{	x = rcauchy(n)
x = cumsum(x)/1:n
plot(1:n, x, type = "l")
abline(h = 0, col = "red")
}


## central limit theorem
central.cauchy = function()
{
  win.graph()
  par(mfrow = c(2,2))
  nt = c(5, 10, 20, 50)
  z = seq(-4, 4, length=100)
  dz = dnorm(z)
  x = xbar = rep(0, 100)
  for (i in 1:4) {
    for (j in 1:100) {
      xbar[j] = sum(rcauchy(nt[i])/nt[i])
      x[j] = sqrt(nt[i])*xbar[j]*(xbar[j]-1)
    }
    plot(density(x), main = paste("sample size =", nt[i]),
         xlim = c(-4,4), ylim = c(0, .5), type = "l")
    par(new = T)
    plot(z, dz, main="", xlab = "", ylab = "", xlim = c(-4,4),
         ylim = c(0, .5), type = "l", col = "red")
    
  }
}
central.cauchy()

## confidence interval
confidence.normal = function(n, nt, mu, std)
{
  win.graph()
  par(pin = c(6,3))
  trial = rep(1, nt)	# trial = I(mu in CI)
  
  t.val = qt(0.975, n-1)
  x = matrix(rnorm(n*nt, mu, std), nrow = nt)
  
  xbar = apply(x, 1, mean)
  xvar = apply(x, 1, var)
  
  limit = t.val * sqrt(xvar/n)
  xbar = xbar - mu
  trial[abs(xbar) > limit] = 0
  trial = cumsum(trial)/1:nt
  plot(1:nt, trial, ylim = c(.9, 1), type = "l")
  abline(h = 0.95, col = "red")
}

confidence.normal(100, 5000, 0, 1)


confidence.beta = function(n, nt, a, b)
{
  win.graph()
  par(pin = c(6,3))
  mu = a/(a+b)
  trial = rep(1, nt)	# trial = I(mu in CI)
  
  t.val = qt(0.975, n-1)
  x = matrix(rbeta(n*nt, a, b), nrow = nt)
  
  xbar = apply(x, 1, mean)
  xvar = apply(x, 1, var)
  
  limit = t.val * sqrt(xvar/n)
  xbar = xbar - mu
  trial[abs(xbar) > limit] = 0
  trial = cumsum(trial)/1:nt
  plot(1:nt, trial, ylim = c(.9, 1), type = "l")
  abline(h = 0.95, col = "red")
}

confidence.beta(100, 5000, 10, 2)
confidence.beta(100, 5000, 6, 6)
confidence.beta(100, 5000, 2, 10)


# regression
reg.confidence = function(nt, n, beta)
{
  win.graph()
  par(mfrow = c(2,1))
  mtitle = c("Intercept", "Slope")
  x = 1:n
  true.y = beta[1] + beta[2] * x
  t.val = qt(0.975, n-2)
  trial.a = trial.b = rep(0, nt)
  
  for (i in 1:nt) {
    y = true.y + rnorm(n)
    lm.out = lm(y ~ x)
    std = sqrt(diag(vcov(lm.out)))
    if (abs(lm.out$coefficients[1]-beta[1]) < t.val*std[1]) trial.a[i]=1
    if (abs(lm.out$coefficients[2]-beta[2]) < t.val*std[2]) trial.b[i]=1
  }
  trial.a = cumsum(trial.a) / 1:nt
  trial.b = cumsum(trial.b) / 1:nt
  plot(trial.a, type = "l", main = mtitle[1], ylim = c(0.6,1))
  abline(0.95,0, col = "red")
  plot(trial.b, type = "l", main = mtitle[2], ylim = c(0.6,1))
  abline(0.95,0, col = "red")
}

beta = c(2, 5)
reg.confidence(1000, 20, beta)


## interval estimation for sd
dice = c(1,2,3,2,6,6,5,1,1,1,4,2,4,1,4,5,6,6,3,2,
         5,6,4,1,2,3,2,2,5,3,5,6,1,4,4,4,3,5,5,1,
         6,1,3,3,2,5,2,2,1,4)

sd(dice)
boot.sample = matrix(0, 200, 50)
std = rep(0, 250)
set.seed(12345)

for(n in 1:200) {
  boot.sample[n,] = sample(dice, replace=T, 50)
  std[n] = sd(boot.sample[n,])
}

summary(std)
hist(std, xlim=c(1.2,2.2), ylim=c(0,50), xlab="bootst.sd")
quantile(std, probs=c(0.025, 0.975))

## jackknife for patch data
data(patch, package="bootstrap")
n = nrow(patch)
y = patch$y
z = patch$z
theta.hat = mean(y)/mean(z)
theta.hat

theta.jack = numeric(n)
for (i in 1:n)
  theta.jack[i] = mean(y[-i])/mean(z[-i])
bias = (n-1) * (mean(theta.jack) - theta.hat)
bias

### EM algorithm

## Example 1: multinomial 
eps = 1e-5
theta = diff = 0.5
k = 0
result = c(k, theta, diff)
while(diff > eps) {
  y = 125 * theta / (theta+2)
  th_hat = (34+y) / (38+34+y)
  diff = abs(theta - th_hat)
  theta = th_hat
  k = k + 1
  result = rbind(result, c(k, theta, diff))
}
round(result, 8)


#지수분포의 log-likelihood 함수값을 구하는 함수
Log.lik = function(x, R=2, lambda, prior)
{   
  lik = 0
  for (r in 1:R)
    lik = lik + prior[r] * dexp(x, rate = lambda[r]) 
  return(sum(log(lik)))
}
#EM 알고리즘
Exp.Mixture = function(X, R=2, maxiter=1000, eps=1e-5)
{
  
  X = as.vector(X)
  N = length(X)
  
  lambda = prior = rep(0, R)
  
  gama = matrix(0, R, N)
  # 확률 초기화
  prior = rep(1/R, R)
  lambda = c(1,3)
  
  old.lik = Log.lik(X, R, lambda, prior) 
  track.lik = as.vector(NULL)
  track.lik = c(old.lik)
  for (i in 1:maxiter)
  {
    for (r in 1:R)
      gama[r, ] = prior[r] * dexp(X, rate = lambda[r])
    denom = apply(gama, 2, sum)
    for (r in 1:R)
    {
      gama[r, ] = gama[r, ] / denom
      lambda[r] = t(gama[r, ]) %*% sqrt(1/X) /  sum(gama[r, ])
    }
    prior = apply(gama, 1, sum) / N
    lambda[1] = 1
    new.lik = Log.lik(X, R, lambda, prior)
    
    if (abs(old.lik - new.lik) < eps * abs(old.lik))  break
    old.lik = new.lik
    track.lik = c(track.lik, old.lik)
  }
  return(list(lambda = lambda, prior = prior, 
              track = track.lik, resp = gama[r, ]))
}

#혼합분포에서 난수 생성
random.number = c()
while(length(random.number)< 100){
    v1 = runif(1)
    if(v1 < 0.1){
      rn = rexp(1, rate = 1)
      random.number = c(random.number, rn)
    }else{
      rn = rexp(1, rate = 5)
      random.number = c(random.number, rn)
    }
}
Exp.Mixture(random.number)


Log.lik = function(x, R=2, mu, sigma, prior)
{   
  lik = 0
  for (r in 1:R)
    lik = lik + prior[r] * dnorm(x, mean = mu[r], sd = sigma[r]) 
  return(sum(log(lik)))
}

Normal.Mixture = function(X, R=2, maxiter=100, eps=1e-5)
{
  X = as.vector(X)
  N = length(X)
  mu = sigma = prior = rep(0, R)
  gama = matrix(0, R, N)
  # find initial centroids using K-means clustering
  prior = rep(1/R, R)
  kmfit = kmeans(X, R)
  mu = kmfit$centers
  sigma = sqrt(kmfit$withinss /(kmfit$size - 1))     
  old.lik = Log.lik(X, R, mu, sigma, prior) 
  track.lik = as.vector(NULL)
  track.lik = c(old.lik)
  for (i in 1:maxiter)
  {
    for (r in 1:R)
      gama[r, ] = prior[r] * dnorm(X, mean = mu[r], sd = sigma[r])
    denom = apply(gama, 2, sum)
    for (r in 1:R)
    {
      gama[r, ] = gama[r, ] / denom
      mu[r] = t(gama[r, ]) %*% X / sum(gama[r, ])
      sigma[r] = sqrt(t(gama[r, ]) %*% (X - mu[r])^2 / sum(gama[r, ]))
    }
    prior = apply(gama, 1, sum) / N
    new.lik = Log.lik(X, R, mu, sigma, prior)
    if (abs(old.lik - new.lik) < eps * abs(old.lik))  break
    old.lik = new.lik
    track.lik = c(track.lik, old.lik)
  }
  return(list(mu = mu, sigma = sigma, prior = prior, 
              track = track.lik, resp = gama[r, ]))
}

Mixture.prob = function(x, mu, sigma, prior)
{
  R = length(mu)
  prob = 0
  for (r in 1:R)
    prob = prob + prior[r] * dnorm(x, mu[r], sigma[r]) 
  return(prob)
}

X = c(-0.39, 0.12, 0.94, 1.67, 1.76, 2.44, 3.72, 4.28, 4.92, 5.53,
      0.06, 0.48, 1.01, 1.68, 1.80, 3.25, 4.12, 4.60, 5.28, 6.22)
fit = Normal.Mixture(X)     

hist(X, nclass = 20, xlab = "X", freq = FALSE, col = 2, ylim = c(0, 1))

X.grid = seq(min(X), max(X), length=100)
plot(X.grid, Mixture.prob(X.grid, fit$mu, fit$sigma, fit$prior), 
     type = "l", ylab = "density", xlab = "X",
     ylim = c(0, 1), col = 2)
lines(sort(X), 1-fit$resp[sort(X, index.return = T, decreasing=T)$ix], 
      type = "b", col = 3)

plot(1:length(fit$track), fit$track, 
     xlab = "iteration", ylab = "Obs. log-likelihood", type = "b", 
     col = 3)



?rexp
n = 100
pi = c(.5, .5)

for(j in 1:n){
  f1 = pi[1] * dexp(random.number, )
}