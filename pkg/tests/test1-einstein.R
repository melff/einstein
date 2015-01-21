library(einstein)


x <- array(seq(3*4*5),c(3,4,5))

x%[%'_i^j_k'

a <- numeric()

a%[%'_i_j_k' <- x%[%'_i^j_k'
str(a)

a <- array(0,c(5,2,2))
a%[%'_k_i_i' <- 1
a[1,,]


m <- matrix(1:9,3,3)
n <- numeric()
m%[%'_i_j' + m%[%'_j_i'

n%[%'_i_j' <- m%[%'_i_j' + m%[%'_j_i'
n

n <- matrix(0,3,3)
n%[%'_i_j' <- m%[%'_i_j' + m%[%'_j_i'
n

n%[%'_i_i' <- n%[%'_i_i' + 2
n

try(m%[%'_i_j' %t% m%[%'^j_i') 

(m%[%'_i_j') %t% (m%[%'^j_i') 

tprod(m%[%'_i_j',m%[%'^j_i') 

