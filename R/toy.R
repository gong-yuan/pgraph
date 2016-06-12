K = 3
n = 100
b = matrix(rnorm(K*n),n,K)
bx = 1:3
by = c(1,2,2)
x = b%*%bx+rnorm(n)
y = b%*%by+rnorm(n)


