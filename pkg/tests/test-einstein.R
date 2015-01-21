library(einstein)


A <- array(seq(3*4*3),c(3,4,3))

array.diag(A,diag.dims=list(c(1,3)))


A <- array(seq(3*4*2*3*4),c(3,4,2,3,4))

# debug(array.diag)
A1 <- array.diag(A,diag.dims=list(c(1,4),c(2,5)))
dim(A1)


x <- tensor(seq(3*4*5),c(3,4,5))

x['_i^j_k']


xx <-  x['_i^j_k']
attributes(xx)

try(xxx <-  xx['_i^j_k'])

u <- tensor()
u['_i^j_k']

arr2mat(x,rows=c(3,2))


xx <- x['_i^j_k'] + x['_i^j_k']
sum(xx)
as.array(x)+as.array(x)

try(xx <- x['_i^j_k'] + x['_i^j_s'])
try(xx <- x['_i^i_k'] + x['_i^i_k'])
try(xx <- x['_i^j_k'] + x['_j^i_k'])

m <- tensor(matrix(1:9,3,3))
m['_i_j'] + m['_j_i']

m['_i_j'] * m['_j_i']
m['_i_j'] %t% m['^j_i'] 
m['_i_j'] %t% m['^j^i'] 
sum(m['_i_j'] * m['^j^i'])
sum(m['_i_j'] * m['^j_i'])

mm <- m['_i_j'] * m['^j_i']

v <- tensor(1:4)

v['_i'] * v['_i']
v['_i'] * v['_j']


m['_i^j'] * m['_j^k']
sum(m['_i^j'] * m['_j^k'])

m1 <- tensor(matrix(1:15,5,3))
m1['_i_j'] * x['^j^k^i']
sum(m1['_i_j'] * x['^j^k^i'])

v['_i'] * v['^i']
sum(v['_i'] * v['^i'])

5*v['^i']

v['_i'] * v['^i'] * v['_i']
sum(v['_i'] * v['^i'] * v['_i'])

m1['_i_j']* x['^j^k^i']*v['_k']

(m1['_i_j']* sum(x['^j^k^i']*v['_k']))
sum(m1['_i_j']* sum(x['^j^k^i']*v['_k']))

sum(m1['_i_j']* x['^j^k^i'])*v['_k']
sum(sum(m1['_i_j']* x['^j^k^i'])*v['_k'])


try(sum.index.permutations(m1['_i_j']))
sum.index.permutations(m['_i_j'])

y <- tensor(array(1,c(3,3,3)))
sum.index.permutations(y['_i_j_k'])
sum.index.permutations(y)

sum.index.partitions(y['_i_j_k'],"i", c("j","k"))


y <- tensor(array(rep(1:3,9),c(3,3,3)))
y['^i_i_k']
sum(y['^i_i_k'])

g <- tensor()

g['_i_j_k']<- y['_j_i_k']


# debug(fill_tensor)
g <- tensor(0,dim(y))
g['_i_j_k']<- y['_j_i_k']
g

g['_i_j_k'] <- m['_j_i']
g
g['_r_s_t']
 
 
g <- tensor(0,c(3,5,3))
g['_i_k_j']<- m['_j_i']
g['_r_s_t']

g <- tensor(0,c(3,5,3))
g['_i_k_i']<- m['_i_i']
g['_r_s_t']

g <- g['_i_j_k']
indices(g)
g['_k_j_i',iperm=TRUE]

a <- array(0,c(3,3,3))
a <- set.array.diag(a,diag.dims=list(c(1,2,3)),1)
a

m["_i_k"]*m["^k_i"]
sum(m["_i_k"]*m["^k_i"])

tprod(m["_i_k"],m["^k_i"])
m["_i_k"] %t% m["^k_i"]
m["_i_k"] %t% m["^k_j"]
m %*% m

v['_i'] * v['^i'] * v['_i']
sum(v['_i'] * v['^i'] * v['_i'])

try(v['_i'] %t% v['^i'] %t% v['_i'])

tensor(a,i="_i_i^k")


a <- array(0,c(3,3,3))
tensor(a)["_i_i^k"] <- 1
a

b <- numeric()
tensor(b)["_i_j_k"] <- tensor(a)["_i_j_k"]
b

a <- tensor(1:9,c(3,3))
o <- tensor(0,c(3,3))
b <- tensor(0:6,c(7,3))
c <- tensor(1:3,3)
a["_i_d"] - b["_t_d"]
b["_t_d"]*c["_d"]
o["_i_d"] - b["_t_d"]*c["_d"]

sum(o["_i_d"] - b["_t_d"]*c["_d"])

sum(o["_i_d"] - b["_t_d"]*c["_d"],over="d")

sum(a["_i_d"],over="d")

diffs <- o["_i_d"] - b["_t_d"]*c["_d"]
sum(diffs^2,over="d")
sum(diffs^2,over=c("d","t"))

diffs <- tensor()
diffs["_t_i_d"] <- o["_i_d"] - b["_t_d"]*c["_d"]
sum(diffs^2,over="d")


X <- array(1,c(6,3,3))
j <- rep(1:2,each=3)
rowsum(X,j)
