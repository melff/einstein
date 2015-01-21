index.partitions <- function(ixlist,show.names=FALSE,verbose=TRUE){
  
  base.set <- unlist(ixlist)
  
  sizes <- sapply(ixlist,length)
  nsets <- length(ixlist)

  ii <- order(sizes,decreasing=TRUE)
  jj <- seq.int(nsets)[ii]
  
  partitions <- listParts(sizes)
  if(verbose) message("created ",length(partitions)," partitions")
  
  partitions <- lapply(partitions,function(x){
    y <- x[jj]
    if(show.names)
      y <- lapply(y,function(y)base.set[y])
    names(y) <- 1:nsets
    structure(y,
              class=c("list","equivalence"))
    })
  partitions
}

# index.partitions(list(
#   c("i"),
#   c("j","k")
# ))
# index.partitions(list(
#   c("i"),
#   c("j","k")
# ),
# show.names=TRUE)
# 
part2perm <- function(x) sapply(x,unlist)
# 
# part2perm(index.partitions(list(
#   c("i"),
#   c("j","k")
# )))

sum.index.partitions <- function(tns,...,verbose=TRUE){
  
  ixlist <- list(...)
  
  d <- dim(tns)
  rnk <- length(d)
  if(length(unique(d))>1) stop("tensor not symmetric")
  
  parts <- index.partitions(ixlist,verbose=verbose)
  prms <- part2perm(parts)
  res <- array(0,dim=d)
  for(i in 1:ncol(prms))
    res <- res + aperm(tns,prms[,i])
  
  structure(
    res,
    indices=attr(tns,"indices"),
    class="tensor")
}

sum.index.permutations <- function(tns){
  
  d <- dim(tns)
  rnk <- length(d)
  if(length(unique(d))>1) stop("tensor not symmetric")
  
  prms <- perms(rnk)
  res <- array(0,dim=dim(tns))
  for(i in 1:ncol(prms))
    res <- res + aperm(tns,prms[,i])
  
  structure(
    res,
    indices=attr(tns,"indices"),
    class="tensor")
}
