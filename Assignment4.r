# Q1
#랜덤 난수 생성
random.number = runif(1000)
#콜모고로프-스미르노프 검정
ks.test(random.number,punif, 0, 1)
#산점도 그리기
plot(random.number)
run.list = c()
for(i in 1:length(random.number) - 1){
  #값이 증가하면 1, 그렇지 않으면 0 추가
  if (isTRUE(random.number[i] < random.number[i+1])) {run.list = append(run.list,1)}
  else {run.list = append(run.list,0)}
}
library(tseries)
#tseries 패키지의 run test
runs.test(factor(run.list))



random.list = random.number
r = 1
for(i in 1:length(run.list) -1){
  if(!isTRUE(run.list[i] == run.list[i+1])){
    r  = r+1 }
}



N = 1000
(r - ((2*N - 1) / 3)) / (sqrt(16*N - 29) / 90)
install.packages('tseries')

length(run.list)
length(random.list)
?sample

# Q2
f = function(x){
  return(cos(x))
}
#적분값의 참값 I
I = integrate(f,0,1)
I
#표본 수는 1000개
N = 1000
hit.or.miss = function(n){
  #column의 개수가 2, row의 개수가 2n개인 행렬 생성
  set.seed(777)
  temp = matrix(runif(2*n), ncol = 2)
  #1번째 Column은 x좌표, 2번때 Column은 Y좌표
  colnames(temp) = c('x', 'y')
  #주어진 영역 내에 점이 들어가는지 확인
  mean(f(temp[,1]) > temp[,2])
}

sample.mean.monte.carlo = function(n){
  set.seed(777)
  #적률법
  random.number = runif(n)
  return(sum(f(random.number)) / n)
}


#중요 함수 생성 - cos(x)의 테일러 전개의 근사식
importance.function = function(x){
  return(6*(1 - 0.5 * x^2) / 5)
}

importance.sampling = function(n){
  set.seed(777)
  random.number = runif(n)
  return(sum(f(random.number) / importance.function(random.number))/n)
}

hit.or.miss(1000)
sample.mean.monte.carlo(1000)
importance.sampling(1000)

g.square = function(x){
  return(cos(x) * cos(x))
}
g.square.div.f = function(x){
  return((cos(x) * cos(x))/importance.function(x))
}
#적중법의 최소 표본 수
2.58 * 2.58 * I$value * (1 - I$value) / 10^(-6)
#표본 평균 몬테칼로 기법의 분산
var.sample.mean = integrate(g.square,0,1)
#표본 평균 몬테칼로 기법의 최소 표본 수
2.58 * 2.58 * (var.sample.mean$value - I$value^2) / 10^(-6)
#주표본기법의 분산
var.importance = integrate(fff, 0, 1)
#주표본기법의 최소 표본 수
2.58 *  2.58 * (var.importance$value - I$value^2) / 10^(-6)

#각 기법에서 최소 표본 수 검산
#적중법
2.58 * sqrt(I$value * (1 - I$value)) / sqrt(887947)
2.58 * sqrt(I$value * (1 - I$value)) / sqrt(887948)
#표본평균 몬테칼로 기법
2.58 * sqrt(var.sample.mean$value - I$value^2) / sqrt(128141)
2.58 * sqrt(var.sample.mean$value - I$value^2) / sqrt(128142)
#주표본 기법
2.58 * sqrt(var.importance$value - I$value^2) / sqrt(1297)
2.58 * sqrt(var.importance$value - I$value^2) / sqrt(1298)



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