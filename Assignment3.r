#Q1

A = matrix(c(1,0,2,-1), nrow = 2)
B = matrix(c(2,-1,1,8,7,-3,6,4), nrow = 2)
x = c(1,3,5,-1)

#(a)
(A + t(A)) %*% B
A %*% B + t(A) %*% B
#(b)
t(B) %*% B
sum(diag(t(B) %*% B))
B %*% t(B)
sum(diag(B %*% t(B)))
#(c)
t(x) %*% t(B) %*% B %*% x
t(B) %*% B %*% x %*% t(x)
sum(diag(t(B) %*% B %*% x %*% t(x)))


#Q2
#힐버트 행렬을 만들어주는 함수
make.Hibert = function(n)
{
  A = matrix(nrow =  n, ncol =  n)
  for(i in 1:n){
    for(j in 1:n){
      A[i,j] = 1/(i+j-1)
    }
  }
  return(A)
}


make.Hibert(2)
det(make.Hibert(2))
solve(make.Hibert(2))


make.Hibert(4)
det(make.Hibert(4))
solve(make.Hibert(4))


make.Hibert(6)
det(make.Hibert(6))
solve(make.Hibert(6))

make.Hibert(8)
det(make.Hibert(8))
solve(make.Hibert(8))

make.Hibert(10)
det(make.Hibert(10))
solve(make.Hibert(10))





#Q3
A = matrix(nrow =  6, ncol =  6)
for(i in 1:6){
  for(j in 1:6){
    A[i,j] = (10 + i - 1)^(j-1)
  }
}
A
B = c(25,16,26,19,21,20)
#역행렬 구하기
aa = solve(A) %*% B

A %*% aa

#LU 분해 - 강의 예제 코드
lufactorization = function(A){
  n = nrow(A)
  L = matrix(0,nrow=n,ncol=n)
  for (k in (1:(n-1))){
    for (i in ((k+1):n)){
      L[i,k] = A[i,k]/A[k,k]
      A[i,k] = 0
      for (j in ((k+1):n)){
        A[i,j] = A[i,j] - L[i,k]*A[k,j]
      }
    }
  }
  for (k in (1:n)){
    L[k,k] = 1
  }
  return(cbind(L,A))
}

lufactorization(A)
aa = solve(lufactorization(A)[,1:6]) %*% B
solve(lufactorization(A)[,7:12]) %*% aa



backsolve(lufactorization(A)[,7:12],forwardsolve(lufactorization(A)[,1:6], B))

y = forwardsolve(t(H3.chol), b) 	# forward elimination (lower)
backsolve(H3.chol,y)			# back substitution (upper)





make.Hibert = function(n)
{
  A = matrix(nrow =  n, ncol =  n)
  for(i in 1:n){
    for(j in 1:n){
      A[i,j] = 1/(i+j-1)
    }
  }
  return(A)
}




#Q4
X = matrix(c(1,2,-4,2,4,2,3,2,1), nrow = 3)
Y = c(1,2,3)
gaussianeliminationpartial(cbind(X, Y))
gauss.X = gaussianeliminationpartial(cbind(X, Y))[,1:3]
gauss.Y = gaussianeliminationpartial(cbind(X, Y))[,4]
solve(gauss.X) %*% gauss.Y
solve(X) %*% Y


gaussianeliminationpartial = function(Ab){
  n = nrow(Ab)
  for (k in (1:(n-1))){
    pivotindex = k
    for (i in ((k+1):n)){
      if (abs(Ab[i,k]) > abs(Ab[pivotindex,k])){
        pivotindex = i
      }
    }
    if (pivotindex != k){
      for (j in (k:(n+1))){
        buffer = Ab[k,j]
        Ab[k,j] = Ab[pivotindex,j]
        Ab[pivotindex,j] = buffer
      }
    }
    for (i in ((k+1):n)){
      mik = Ab[i,k]/Ab[k,k]
      Ab[i,k] = 0
      for (j in ((k+1):(n+1))){
        Ab[i,j] = Ab[i,j] - mik*Ab[k,j]
      }
    }
  }
  return(Ab)
}

H3.chol = chol(H3)
H3.chol
crossprod(H3.chol, H3.chol)
chol2inv(H3.chol)	# inverse

#Q5 - (a)
# Choleski decomposition
X = matrix(c(1,1,1,1,1,1,1,2,3,5,5,7,1,3,3,4,4,5), nrow = 6)
Y = c(2,4,5,8,8,9)
choleskyfactorization(t(X) %*% X)
t(choleskyfactorization(t(X) %*% X))
chol(t(X) %*% X)

lm(Y ~ X)$coef

#강의 예제 함수
choleskyfactorization = function(A){
  n = nrow(A)
  L = matrix(0,nrow=n,ncol=n)
  for (i in (1:n)){
    L[i,i] = A[i,i]
    if (i > 1){
      for (k in (1:(i-1))){
        L[i,i] = L[i,i] - L[i,k]*L[i,k]
      }
    }
    L[i,i] = (L[i,i])^(1/2)
    if (i < n){
      for (j in ((i+1):n)){
        L[j,i] = A[j,i]
        if (i > 1){
          for (k in (1:(i-1))){
            L[j,i] = L[j,i] - L[j,k]*L[i,k]
          }
        }
        L[j,i] = L[j,i]/L[i,i]
      }
    }
  }
  return(L)
}

#Q5 - (b) : Regression with QR Decomposition
X = matrix(c(1,1,1,1,1,1,1,2,3,5,5,7,1,3,3,4,4,5), nrow = 6)
Y = c(2,4,5,8,8,9)
X.qr = qr(X)
qr.Q(X.qr)
qr.R(X.qr)
solve(qr.R(X.qr)) %*% t(qr.Q(X.qr)) %*% Y
qr.solve(X, Y)
lm(Y ~ X)$coef



# QR decomposition
X.qr = qr(X)
diag(X.qr$qraux, 3, 3)
Q = qr.Q(X.qr)
Q
R = qr.R(X.qr)
R
Q %*% R
qr.solve(R) %*% t(Q)	# inverse


#Q6
A = matrix(c(2,-4,2,4,-9,7,2,1,3), nrow = 3, byrow = T)
t(A) %*% A
eigen(t(A) %*% A, symmetric = T)
svd(A)
EA = eigen(A, symmetric = T)
(svd(A)$d)^2
det(A)
prod(EA$values)
# spectral decomposition
A = matrix(c(5,25,35,25,155,175,35,175,325), ncol = 3)
A
EA = eigen(A, symmetric = T)
EA
det(A)
prod(EA$values)
H3.svd = svd(H3)
H3.svd
H3.svd$u %*% diag(H3.svd$d) %*% t(H3.svd$v)
H3.svd$v %*% diag(1/H3.svd$d) %*% t(H3.svd$u) # inverse




# constructing matrix
H3 = 1/cbind(1:3, 2:4, 3:5)
matrix(seq(1,12), nrow = 3)
x = seq(1,3)
x2 = x^2
X = cbind(x, x2)
X
X = matrix(c(1,2,3,1,4,9), ncol = 2)
X

# accessing matrix elements
X[3,2]
X[3,]
X[,2]
colnames(X)
rownames(X)
rownames(X) = c("obs1","obs2","obs3")
X

# matrix properties
mat = rep(1:4, rep(3,4))
dim(mat) = c(3,4)
mat
dim(mat)

# diagonal matrix
D = diag(c(1,2,3))
D
diag(D)
I = diag(3)
I

# triangluar matrix
lower.tri(H3)
Hnew = H3
Hnew[upper.tri(H3)] = 0
Hnew

# transpose
A = matrix(c(rep(1,5), -1, 1, -1), byrow = T, ncol = 4)
A
A1 = t(A)
A1

# matrix arithmetic
Y = 2 * X
Y
Y + X
t(Y) + X
X * Y

# matrix multiplication
t(Y) %*% X
Y %*% X
crossprod(Y,X)

# determinant and inverse
det(H3)
H3inv = solve(H3)
H3inv
H3inv %*% H3 

# generalized inverse
library(limSolve)
A = matrix(c(1:8, 6, 8, 10, 12), nrow=4, ncol=3)
B = 0:3
X = Solve(A, B)
A %*% X - B
(gA = Solve(A))



# solving linear equations
A = matrix(c(4,2,1,-2,-1,3,3,1,-1), ncol = 3)
A
b = c(9,3,4)
solve(A,b)

# Singular value decomposition
H3.svd = svd(H3)
H3.svd
H3.svd$u %*% diag(H3.svd$d) %*% t(H3.svd$v)
H3.svd$v %*% diag(1/H3.svd$d) %*% t(H3.svd$u) # inverse

# Choleski decomposition
H3.chol = chol(H3)
H3.chol
crossprod(H3.chol, H3.chol)
chol2inv(H3.chol)	# inverse
# H_3 x = b, b = (1,2,3)'
b = 1:3
y = forwardsolve(t(H3.chol), b) 	# forward elimination (lower)
backsolve(H3.chol,y)			# back substitution (upper)

# QR decomposition
H3.qr = qr(H3)
H3.qr
Q = qr.Q(H3.qr)
Q
R = qr.R(H3.qr)
R
Q %*% R
qr.solve(R) %*% t(Q)	# inverse

# LU decomposition
set.seed(1)
library(Matrix)
mm = Matrix(round(rnorm(9),2), nrow=3)
mm

lum = lu(mm)
str(lum)

elu = expand(lum)
elu

# spectral decomposition
A = matrix(c(5,25,35,25,155,175,35,175,325), ncol = 3)
A
EA = eigen(A, symmetric = T)
EA
det(A)
prod(EA$values)


# outer product
x1 = 1:5
outer(x1, x1, "/")
outer(x1, x1, "-")
y = 5:10
outer(x1, y, "+")



lufactorization = function(A){
  n = nrow(A)
  L = matrix(0,nrow=n,ncol=n)
  for (k in (1:(n-1))){
    for (i in ((k+1):n)){
      L[i,k] = A[i,k]/A[k,k]
      A[i,k] = 0
      for (j in ((k+1):n)){
        A[i,j] = A[i,j] - L[i,k]*A[k,j]
      }
    }
  }
  for (k in (1:n)){
    L[k,k] = 1
  }
  return(cbind(L,A))
}
diag(c(1,2,3),3,3)
