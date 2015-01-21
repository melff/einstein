
collapse.array <- function(x,dims=1){
  
  if(!length(dims)) return(x)
  sdx <- seq_along(dim(x))
  ksdx <- sdx[-dims]
  dsdx <- sdx[dims]
  if(length(ksdx))
    rowSums(aperm(x,c(ksdx,dsdx)),
            dims=length(ksdx))
  else sum(x)
}


check.eq.dim <- function(tocheck,d) length(unique(d[tocheck]))>0
elems <- function(i,x)x[i]

array.diag <- function(x,diag.dims=NULL){
  
  if(is.null(diag.dims)) return(x)
  d <- dim(x)
    
  diag.dims.ok <- sapply(diag.dims,check.eq.dim,d=d)
  if(!all(diag.dims.ok)) stop("cannot extract diagonal of non-equal dimensions")
  
  sizes <- lapply(diag.dims,elems,d)
  ds <- sapply(sizes,"[",1)
  
  orders <- lapply(diag.dims,length)

  krons <- mapply(kron,orders,ds,SIMPLIFY=FALSE)
  d.use <- lapply(krons,as.logical)
  
  dd <- unlist(diag.dims)
  rem <- seq_along(d)
  rem <- rem[!(rem%in%dd)]
  
  perm <- c(rem,dd)
  y <- aperm(x,perm)
  if(length(rem)) pdrem <- prod(d[rem])
  else pdrem <- 1
  
  psz <- sapply(sizes,prod)
  dy <- c(pdrem,psz)
  
  dim(y) <- dy
  
  ii <- rep(list(TRUE),length(rem))
  ii <- c(ii,d.use)
  
  do.call("[",c(list(y),ii))
}

set.array.diag <- function(x,diag.dims,value){
  if(is.null(diag.dims)) return(x)
  d <- dim(x)
  
  diag.dims.ok <- sapply(diag.dims,check.eq.dim,d=d)
  if(!all(diag.dims.ok)) stop("cannot assign to diagonal of non-equal dimensions")
  
  sizes <- lapply(diag.dims,elems,d)
  ds <- sapply(sizes,"[",1)
  
  orders <- lapply(diag.dims,length)
  
  krons <- mapply(kron,orders,ds,SIMPLIFY=FALSE)
  d.use <- lapply(krons,as.logical)
  
  dd <- unlist(diag.dims)
  rem <- seq_along(d)
  rem <- rem[!(rem%in%dd)]
  
  rem.dim <- dim(x)[rem]
  
  perm <- c(dd,rem)
  y <- aperm(x,perm)
  
  psz <- sapply(sizes,prod)
  dim(y) <- c(psz,rem.dim)
  
  ii <- rep(list(TRUE),length(rem))
  ii <- c(d.use,ii)
  
  y <- do.call("[<-",c(list(y),ii,list(value=value)))
  dim(y) <- c(unlist(sizes),rem.dim)
  
  rev.perm <- seq_along(perm)
  rev.perm[perm] <- rev.perm
  
  aperm(y,rev.perm)
}

kron <- function(order=2,size){
  
  stopifnot(length(size)==1)
  if(order==2)return(diag(nrow=size))
  dm <- rep(size,order)
  k <- array(0,dim=dm)
  i <- as.matrix(1:size)
  iii <- i[,rep(1,order)]
  k[iii] <- 1
  k
}

n.all.eq <- function(x) length(unique(x))>1




grid.matrix <- function(...){
  
  arglist <- c(...)
  
  res <- numeric(0)
  
  for(arg in arglist){
    if(NROW(res)==0) res  <- as.matrix(seq.int(arg))
    else {
      m <- nrow(res)
      i <- rep(1:m,arg)
      res <- cbind(res[i,],rep(1:arg,each=m))
    }
  }
  res
}

index.array <- function(a) grid.matrix(dim(a))

rowsum.array <- function(x, group, reorder=TRUE, ...){
  
  dx <- dim(x)
  dim(x) <- c(dx[1],prod(dx[-1]))
  x <- rowsum(x,group,reorder,...)
  dim(x) <- c(dim(x)[1],dx[-1])
  x
}

