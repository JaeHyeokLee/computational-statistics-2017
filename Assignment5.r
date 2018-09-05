#Q1. Pareto Distribution
u = runif(10000)
Fx = function(x){
  return(1 - (2/x)^2)
}
F.inverse = function(x){
  return(2/sqrt(1-x))
}
F.density = function(x){
  return(8/x^3)
}
hist(F.inverse(u), xlim = c(2, 3), freq = FALSE, breaks = 10000)
curve(F.density, add = TRUE)
library(lattice)
qqmath(F.inverse(u),distribution =  function(p) 2/sqrt(1-p), xlim = c(0, 60), ylim = c(0, 60))
#qpareto함수를 사용하기 위해 rmutil 패키지를 이용했습니다.
#qpareto 함수를 이용해도 같은 결과를 얻을 수 있습니다.
library(rmutil)
qqmath(F.inverse(u),distribution =  function(p) qpareto(p, 2, 2), xlim = c(0, 60), ylim = c(0, 60))

qqplot(quantile(F.inverse(u), seq(0, 1, by = 0.01)), F.density(u))
qqmath~F.inverse(u), distribution = F.density)
?qqmath
??qpareto
??quantile
install.packages('rmutil')

?qpareto
plot(F.inverse(u))
qqmath(F.inverse(u))

??lattice
qqline(F.inverse(u))
quantile(F.inverse(u))[2]
quantile(F.inverse(u))[4]
?qqline
qqmath(quantile(F.inverse(u), seq(0,1,0.001)),distribution =  F.density)
?seq
?qchisq
qq(F.inverse(u) ~ F.density(u))
plot(F.density)
#Q2. Random Number with Logistic
random.number = runif(10000, -1, 0)
logistic = log(-1-1/random.number)
f = function(x){
  return(exp(x) / (1 + exp(x))^2)
}
hist(logistic, xlim = c(-2, 2), freq = FALSE, breaks = 200)
curve(f, add = TRUE)
logistic



#Q4. Box-Muller Transform
box.muller = function(n){
  #n이 짝수인 경우에만 사용 가능  
  u1 = runif(n/2)
  u2 = runif(n/2)
  z1 = sqrt(-2*log(u1)) * cos(2*pi*u2)
  z2 = sqrt(-2*log(u1)) * sin(2*pi*u2)
  return(cbind(z1, z2))
}
#Q4. Polar Method
polar = function(n){
  random.number = c()
  i = 0
  while(length(random.number)< n){
    i = i+1
    v1 = runif(1, -1, 1)
    v2 = runif(1, -1, 1)
    w.square = v1^2 + v2^2
    if(w.square < 1){
      z1 = v1 * sqrt(-1*log(w.square) / w.square)
      z2 = v2 * sqrt(-1*log(w.square) / w.square)  
      random.number = c(random.number, z1, z2)
    }
  }
  print(i)
  return(random.number)
}


system.time(box.muller(100000))
head(box.muller(100000), 10)
system.time(polar(100000))
head(polar(100000), 10)


get.v = function(){
  v1 = runif(1, -1, 1)
  v2 = runif(1, -1, 1)
  if(v1^2 + v2^2 < 1){
    return(cbind(v1, v2))
  }
}

get.v()
gamma(4)
polar = function(n){
  random.number = c()
  i = 0
  while(2 * length(random.number)/3 < n){
      i = i+1
      v1 = runif(1, -1, 1)
      v2 = runif(1, -1, 1)
      w.square = v1^2 + v2^2
      if(w.square < 1){
        random.number = rbind(random.number, cbind(v1, v2, w.square))  
        }
  }
  z1 = sqrt(-2*log(random.number[,3]) / random.number[,3]) * random.number[,1]
  z2 = sqrt(-2*log(random.number[,3]) / random.number[,3]) * random.number[,2]
  print(i)
  return(cbind(z1, z2))
}





polar1 = function(n){
  result = c()
  i = 0
  while(length(result) < n){
    i = i+1
    v1 = runif(1,-1,1)
    v2 = runif(1,-1,1)
    
    w = sqrt(v1^2 + v2^2)
    
    if(w > 1) next
    
    z1 = v1 * sqrt(-1*log(w^2) / w^2)
    z2 = v2 * sqrt(-1*log(w^2) / w^2)
    
    result = c(result, z1, z2)

  }
  print(i)
  return(result)
}
#사용자  시스템 elapsed 
#46.17    0.17   48.33

#사용자  시스템 elapsed 
#45.69    0.19   47.64 

#사용자  시스템 elapsed 
#46.55    0.17   49.61


system.time(polar(100000))

#사용자  시스템 elapsed 
#74.86    0.78   77.50
#63496


#[1] 232621
#사용자  시스템 elapsed 
#62.34    0.51   65.03

#사용자  시스템 elapsed 
#34.03    0.23   35.49 

system.time(polar1(100000))


#63782 w가 큰걸로 했을때는 이정도.
#233864
runif(2, -1, 1)
?runif

aa = polar1(100)
length(aa)
aa[,1]
aa[,2][1]^2
aa[,3][1]
v1 = runif(1, -1, 1)
v2 = runif(1, -1, 1)
a = c()
length(a)
a = rbind(a, cbind(v1, v2))
a = rbind(a, cbind(v1, v2))
a
a[,2]
a[1,]

#Q3. Beta Distribution
#기각법을 이용해서 Gamma Distribution 에서 난수 생성
#교재에 있는 코드를 참고했습니다.
rgamma = function(n, alpha, beta){
  rgamm = rep(0, n)
  if(alpha > 1){
    lambda = sqrt(2*alpha - 1)
    aln4 = alpha - log(4)
    delta = alpha + lambda
    for(i in 1:n){
      repeat{
        uniform = runif(2)
        v = log(uniform[1] / (1 - uniform[1])) / lambda
        rgamm[i] = alpha * exp(v)
        if(log(uniform[1]*uniform[1]*uniform[2]) ==
           (aln4 + delta * v - rgamm[i])) break
      }
    }
  }
  else{
    e = exp(1)
    kappa = (e + alpha) / e
    for(i in 1:n){
      repeat{
        uniform = runif(2)
        p = uniform[1] * kappa
        if(p>1){
          rgamm[i] = -log((kappa - p) / alpha)
          if(uniform[2] < rgamm[i]^(alpha-1)) break
        }
        else{
          rgamm[i] = p^(1/alpha)
          if(uniform[2] <= exp(-x)) break
        }
      }
    }
  }
  return(beta * gamma)
}
#감마분포에서 만든 난수로 베타분포 난수 생성
gamma.to.beta = function(n, alpha, beta){
  x1 = rgamma(n, alpha, 1)
  x2 = rgamma(n, beta, 1)
  return(x1/(x1 + x2))
}
#h(x) = ax^(a - 1) 로 설정한 겅우 기각법으로 난수를 생성하는 함수
beta.rejection = function(n, alpha, beta){
  h = function(x){
    return(alpha * x^(alpha - 1))
  }
  H = function(x){
    return(x^(alpha))
  }
  H.inverse = function(x){
    return(x^(1/alpha))
  }
  random.number = c()
  while(length(random.number) < n){
    uniform = runif(2)
    inverse = H.inverse(uniform[1])
    #g(x) = ((1-x)^(beta - 1)) / alpha
    temp = ((1 - inverse)^(beta - 1))/alpha
    if(uniform[2] < temp){
      random.number = c(random.number, inverse)
    }
  }
  return(random.number)
}

gamma.to.beta(10, 2, 3)
beta.rejection(10, 2, 3)

t1 = rgamma(10000, 2, 1)
t2 = rgamma(10000, 5, 1)
tt = t1/(t1 + t2)

gamma11 = gamma.rejection(10000, 2, 5)
gamma11 = gamma.rejection(10, 2, 5)
ks.test(tt, gamma11)

# pseudo random number
runif(5)
runif(10, min=-3, max=-1)
set.seed(32789)
runif(5)
set.seed(32789)
runif(5)
??random
# Bernoulli r.v.
set.seed(23207)
guesses = runif(20)
correct.answers = (guesses < 0.2)
correct.answers
table(correct.answers)

# Binomial r.v.
dbinom(x=4, size=6, prob=0.5)		# Pr(X=4)
pbinom(4, 6, 0.5)				# Pr(X <= 4)
qbinom(0.89, 6, 0.5)			# 89 percentile

defectives = rbinom(24, 15, 0.1)
defectives
any(defectives > 5)

# Poisson r.v.

dpois(x=3, lambda=0.5)
rpois(10, 3.7)

# exponential r.v.
pexp(1, rate = 3)

# normal r.v.
qnorm(0.95, mean=2.7, sd=3.3)
rnorm(10,-3,0.5)
# x ~ N(0,1) conditional on 0<x<3
x = rnorm(10000)
x = x[(0<x) & (x<3)]
hist(x, probability=T)

# multivariate normal r.v.
library(MASS)
mu = c(0, 1)
sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2, byrow = T)
mvrnorm(10, mu, sigma)


## Hit-or-Miss Monte Carlo
# Example 6-1
myPi.1 = function(n)
{
  x = matrix(runif(n*2), ncol=2)
  r = mean(apply(x^2, 1, sum) < 1)
  return(4*r)
}

myPi.1(4000)
x = matrix(runif(10), ncol=2)
mean(apply(x^2, 1, sum) < 1)
## sample mean Monte Carlo
# pi
Sample.Mean = function(n) 
{
  x = runif(n)
  return(4*mean(sqrt(1-x^2)))
}

Sample.Mean(4000000)



## importance sampling
# int_0^1 exp(-x)/(1+x^2) dx

m = 10000
theta.hat = se = numeric(3)
g = function(x) {
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}

# candidate 1: I(0<x<1)
x = runif(m)     
fg = g(x)
theta.hat[1] = mean(fg)
se[1] = sd(fg)

# candiate 2: exp(-x)
x = rexp(m, 1) 
fg = g(x) / exp(-x)
theta.hat[2] = mean(fg)
se[2] = sd(fg)

# candiate 3: 4*(pi *(1+x^2))^{-1} I(0<x<1) (Cauchy dist on (0,1))
u = runif(m)    # inverse transform method
x = tan(pi * u / 4)
fg = g(x) / (4 / ((1 + x^2) * pi))
theta.hat[3] = mean(fg)
se[3] <- sd(fg)

rbind(theta.hat, se)



# Buffon's needle
Buffon = function(n, lofneedle, distance)
{
  lofneedle = lofneedle / 2
  distance = distance / 2
  r1 = runif(n)
  r2 = runif(n)
  prob = mean(r1*distance < lofneedle*sin(r2*pi))
  return(prob)
}

Buffon(5000,15,20)


U = runif(100000, 1, 7)		# multiple integration
V = runif(100000, 3, 10)
mean(sin(U-V))*42		# 42: density
X = rexp(100000)


## Inverse transform
rexponential = function(n, lambda)
{
  if (lambda <= 0) stop("lambda should be positive")
  return(-log(runif(n))/lambda)
}

rexponential(10, 0.5)

rBernoulli = function(n, p)
{
  if (p<0 || p>1) stop("p should be in [0,1]")
  q = 1 - p
  return(ceiling(runif(n)-q))
}

# Rejection
rnormal = function(n, mu=0, std=1)
{
  if (std<=0) stop("std should be positive")
  r = rep(0, n)
  for (i in 1:n) 
  {
    repeat
    {
      u = runif(2)
      u = - log(u)
      if (u[2] > (u[1]-1)^2)
      {
        r[i] = u[1]
        if (runif(1) < 0.5) r[i] = -r[i]
        break
      }
    }
    r = std*r + mu
  }
  return(r)
}

a = runif(2)  
a[2]
rbinomial.rejection = function(n, size, prob)
{
  if (prob<0 || prob>1) stop("prob must be in [0,1]")
  p = ifelse(prob<=0.5, prob, 1-prob)
  a = sqrt(2*size*p*(1-p))
  x0 = size*p
  rbinom = rep(0, n)
  for (i in 1:n)
  {
    repeat
    {
      u12 = runif(2)
      yx0 = a * tan(pi*u12[1])
      y = yx0 + x0
      gx = 1.2*a*(1+yx0^2/(a^2))*dbinom(as.integer(y), size, p)
      if (u12[2] <= gx) break
    }
    rbinom[i] = as.integer(y)
  }
  rbinom = ifelse(prob<=0.5, rbinom, size - rbinom)
  
  return(rbinom)
}

## discrete distributions
rdiscunif = function(n, a, b)
{
  if (as.integer(a) == n && as.integer(b) == b)
    return(as.integer(a+(b-a+1)*runif(n)))
  else stop("a and b should be integers.")
}

rbinomial = function(n, size, prob)
{
  if (prob>1 || prob<0) stop("p should be in [0,1]")
  p = ifelse(prob<=0.5, prob, 1-prob)
  f0 = (1-p)^size
  p1p = p/(1-p)
  ranbin = runif(n)
  for (i in 1:n)
  {
    x = 0
    fx = f0
    repeat
    {
      if (ranbin[i] <= fx)
      {
        ranbin[i] = x
        break
      } else
      {
        ranbin[i] = ranbin[i] - fx
        fx = (size-x)/(x+1)*p1p*fx
        x = x + 1
      }
    }
  }
  ranbin = ifelse(prob<=0.5, ranbin, size-ranbin)
  return(ranbin)
}

rpoisson = function(n, lambda)
{
  ep = exp(-lambda)
  rpoiss = rep(0, n)
  for (i in 1:n)
  {
    tr = 1
    repeat
    {
      tr = tr*runif(1)
      if (tr <= ep) break
      else rpoiss[i] = rpoiss[i] + 1
    }
  }
  return(rpoiss)
}


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
central.Exponential = function(lambda)
{
  win.graph()
  par(mfrow = c(2,2))
  nt = c(5, 10, 20, 50)
  z = seq(-4, 4, length=100)
  dz = dnorm(z)
  x = xbar = rep(0, 100)
  for (i in 1:4) {
    for (j in 1:100) {
      xbar[j] = sum(rexp(nt[i], lambda))/nt[i]
      x[j] = sqrt(nt[i])*xbar[j]*(xbar[j]-1/lambda)
    }
    plot(density(x), main = paste("sample size =", nt[i]),
         xlim = c(-4,4), ylim = c(0, .5), type = "l")
    par(new = T)
    plot(z, dz, main="", xlab = "", ylab = "", xlim = c(-4,4),
         ylim = c(0, .5), type = "l", col = "red")
    
  }
}

central.Exponential(1)

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