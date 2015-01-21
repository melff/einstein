indices <- function(x)  attr(x,"indices")
index.names <- function(x) {
  ix <- attr(x,"indices")
  substr(ix,2,nchar(ix))
}
index.types <- function(x) {
  ix <- attr(x,"indices")
  substr(ix,1,1)
} 
is.subscript <- function(x) substr(x,1,1)=="_"
is.superscript <- function(x) substr(x,1,1)=="^"

split.indices <- function(i){
  inames <- unlist(strsplit(i,"[_\\^]"))[-1]
  itypes <- regmatches(i,gregexpr("[_\\^]",i))[[1]]
  if(length(inames)!=length(itypes)) 
    stop("malformed index expression")
  list(names=inames,types=itypes)
}

same.type <- function(x,y) substr(x,1,1)==substr(y,1,1)

combine.type <- function(t){
  ut <- unique(t)
  if(length(ut)>1) "+" else ut
}

add_indices <- function(x,i,override=TRUE,iperm=FALSE){
  
  it <- split.indices(i)
  inames <- it$names
  xnames <- index.names(x)
  
  if(iperm){
    if(anyDuplicated(inames))stop("duplicate indices not allowed for permutation")
    if(override)stop("cannot permute indices and override at the same time")
    if(!length(xnames))stop("no indices to permute")
    if(is.array(x)){
      perm <- match(inames,xnames)
      if(any(is.na(perm)))stop("cannot permute for undefined index")
      x <- aperm(x,perm)
    }  
  } else if(!override && length(xnames)) stop("indices already defined")

  ixs <- paste0(it$types,it$names)

  if(length(x)){
    
    if(!is.array(x)) 
      dim(x) <- length(x)
    
    if(length(dim(x)) != length(ixs)) stop("number of indices must match number of dims")
    
    perhaps.diag(structure(x,
                               indices=ixs,
                               class= "tensor"
    ))
  }
  else {
    if(anyDuplicated(substr(ixs,2,nchar(ixs)))) 
      stop("duplicate indices not allowed for for tensor skeleton")
    structure(x,
              indices = ixs,
              class="tensor_skeleton")
  }
}

"[.tensor" <- function(x,i,...,iperm=FALSE,override=!iperm){
  if(is.character(i)) add_indices(x,i,override=override,iperm=iperm)
  else unclass(x)[i,...]
}
tensor <- function(x=numeric(),dim=NULL,indices=NULL){
  
  if(length(dim))
    x <- array(x,dim)
  else if(!is.array(x))
    dim(x) <- length(x)

  if(length(indices))
    add_indices(x,i=indices,override=TRUE)
  else structure(x,indices=NULL,
                 class=if(length(x))"tensor" else "tensor_skeleton")
}

"tensor<-" <- function(x,value){
  if(inherits(x,"tensor"))
    invisible(value)
  else 
    {
      invisible(drop(as.array(value)))
    }
}

"[<-.tensor" <- function(x,i,value) {
  if(inherits(value,"tensor")) 
    x <- fill_tensor_tensor(x,i,value)
  else 
    x <- fill_tensor_numeric(x,i,value)
  invisible(x)
}

"[<-.tensor_skeleton" <- function(x,i,value) {
  x <- fill_tensor_skeleton(x,i,value)
  invisible(x)
}

arr2mat <- function(x,rows,cols){
  
  d <- seq_along(dim(x))
  
  if(missing(rows)&&missing(cols)) {
    rows <- 1
    cols <- setdiff(d,rows)
  }
  else if(missing(rows)) rows <- setdiff(d,cols)
  else if(missing(cols)) cols <- setdiff(d,rows)
  
  aperm.dims <- c(rows,cols)
  if(length(unique(aperm.dims))!=length(dim(x)))
    stop("'rows' and 'cols' must match dim of argument")
  
  matrix.dims <- c(prod(dim(x)[rows]),
                   prod(dim(x)[cols]))
  
  y <- aperm(x,aperm.dims)
  dim(y) <- matrix.dims
  return(y)
}





"Ops.tensor" <- function(e1,e2){
  
  d1 <- dim(e1)
  d2 <- dim(e2)
  
  i1 <- indices(e1)
  i2 <- indices(e2)
  if(!length(i1)) {
    if(inherits(e1,"tensor"))
      stop("left-hand size does not have index names")
    else if(length(e1)>1) stop("left-hand side must be a tensor or a scalar")
    NextMethod(e1,e2)
  }
  if(!length(i2)) {
    if(inherits(e2,"tensor"))
      stop("right-hand size does not have index names")
    else if(length(e2)>1) stop("right-hand side must be a tensor or a scalar")
    NextMethod(e1,e2)
  }
  
  if(anyDuplicated(i1) || anyDuplicated(i2))
    stop("indices must be unique")
  
  if(length(i1)==length(i2) && all(i1 %in% i2)){
    
    #message("elementwise operation")
    if(!all(i1==i2)){
      imatch <- match(i2,i1)
      e2 <- aperm(e2,imatch)
    }
    perhaps.diag(structure(NextMethod(),
              indices=i1,
              class="tensor"))
  }
  else {
  
    in1 <- index.names(e1)
    in2 <- index.names(e2)
    itps1 <- index.types(e1)
    itps2 <- index.types(e2)
    
    common.names <- intersect(in1,in2)
    if(!length(common.names)){
      #message("outer product")
      d <- c(dim(e1),dim(e2))
      ixs <- c(i1,i2)
      res <- outer(as.vector(e1),as.vector(e2),.Generic)
      dim(res) <- d
      perhaps.diag(structure(res,
                indices=ixs,
                class="tensor"))
    }
    else{
      
      ii1 <- match(common.names,in1)
      ii2 <- match(common.names,in2)
      
      uq1 <- which(!(in1 %in% common.names))
      uq2 <- which(!(in2 %in% common.names))
      
      e1 <- arr2mat(e1,cols=uq1,rows=ii1)
      d1u <- d1[uq1]
      d1c <- d1[ii1]
      
      e2 <- arr2mat(e2,cols=uq2,rows=ii2)
      d2u <- d2[uq2]
      d2c <- d2[ii2]
      
      if(!all(d1c==d2c)) stop("dimensions do not match")
      
      e2 <- array(e2,c(dim(e2),ncol(e1)))
      e1 <- as.vector(e1)
      e2 <- aperm(e2,c(1,3,2))
      res <- NextMethod()
      
      dim(res) <- c(d1c,d1u,d2u)

      uq1.nms <- in1[uq1]
      uq2.nms <- in2[uq2]
      r.nms <- c(common.names,uq1.nms,uq2.nms)

      ctps1 <- itps1[ii1]
      ctps2 <- itps2[ii2]
      uq1.tps <- itps1[uq1]
      uq2.tps <- itps2[uq2]
      
      common.types <- ifelse(ctps1==ctps2,ctps1,"+")
      r.types <- c(common.types,uq1.tps,uq2.tps)    
      ixs <- paste0(r.types,r.nms)
      
      in12 <- unique(c(in1,in2))
      perm <- match(in12,r.nms)
      res <- aperm(res,perm)
      ixs <- ixs[perm]
      
      perhaps.diag(structure(res,
                indices=ixs,
                class="tensor"))
      
    }
  }
}



"%t%" <- tprod <- function(e1,e2){
  #message("tprod")
  
  d1 <- dim(e1)
  d2 <- dim(e2)
  
  i1 <- indices(e1)
  i2 <- indices(e2)
  if(!length(i1)) 
      stop("left-hand size does not have index names")
  if(!length(i2)) 
    stop("right-hand size does not have index names")
  
  if(anyDuplicated(i1) || anyDuplicated(i2))
    stop("indices must be unique")
  
  if(length(i1)==length(i2) && all(i1 %in% i2)){
    
    #message("elementwise product")
    if(!all(i1==i2)){
      imatch <- match(i2,i1)
      e2 <- aperm(e2,imatch)
    }
    structure(unclass(e1)*unclass(e2),
              indices=i1,
              class=c("elementwise tensor product","tensor"))
  }
  else {
    
    in1 <- index.names(e1)
    in2 <- index.names(e2)
    itps1 <- index.types(e1)
    itps2 <- index.types(e2)
    
    common.names <- intersect(in1,in2)
    if(!length(common.names)){
      #message("outer product")
      d <- c(dim(e1),dim(e2))
      ixs <- c(i1,i2)
      res <- as.vector(e1)%o%as.vector(e2)
      dim(res) <- d
      structure(res,
                indices=ixs,
                class=c("outer tensor product","tensor"))
    }
    else{
      
      ii1 <- match(common.names,in1)
      ii2 <- match(common.names,in2)
      cix1 <- i1[ii1]
      cix2 <- i2[ii2]
      
      uq1 <- which(!(in1 %in% common.names))
      uq2 <- which(!(in2 %in% common.names))
      
      cix1.cix2.diff.type <- !same.type(cix1,cix2)
      
      if(!any(cix1.cix2.diff.type)){
        
        #message("(partial) outer product")
        
        e1 <- arr2mat(e1,cols=uq1,rows=ii1)
        d1u <- d1[uq1]
        uqix1 <- i1[uq1]
        d1c <- d1[ii1]
        
        e2 <- arr2mat(e2,cols=uq2,rows=ii2)
        d2u <- d2[uq2]
        uqix2 <- i2[uq2]
        d2c <- d2[ii2]
        
        if(!all(d1c==d2c)) stop("dimensions do not match")
        
        res <- array(e2,c(dim(e2),ncol(e1)))
        res <- aperm(res,c(1,3,2)) * as.vector(e1)
        
        dim(res) <- c(d1c,d1u,d2u)
        ixs <- c(cix1,uqix1,uqix2)
        
        in12 <- unique(c(in1,in2))
        perm <- match(in12,r.nms)
        res <- aperm(res,perm)
        ixs <- ixs[perm]
        
        structure(res,
                  indices=ixs,
                  class=c("(partial) outer tensor product","tensor"))
        
      }
      else {
        
        cin1 <- in1[ii1]
        cin2 <- in2[ii2]
        clps <- cin1[cix1.cix2.diff.type]
        
        rem1 <- which(!(in1 %in%clps))
        rem2 <- which(!(in2 %in%clps))
        
        clps1 <- which(in1 %in%clps)
        clps2 <- match(in1[clps1],in2)
        
        e1 <- arr2mat(e1,rows=rem1,cols=clps1)
        e2 <- arr2mat(e2,cols=rem2,rows=clps2)
        res <- drop(e1 %*% e2)
        
        d.res <- c(d1[rem1],d2[rem2])
        
        if(!length(d.res)) return(res)
        dim(res) <- d.res
        ixs <- c(i1[rem1],i2[rem2])
        
        perhaps.diag(structure(res,
                  indices=ixs,
                  class=c("(partial) inner tensor product","tensor")))
        
      }
    }
  }
}

as.array.tensor <- function(x,iperm=NULL,...){

  if(length(iperm)) x <- tperm(x,iperm)
  idxs <- indices(x)
  inms <- index.names(x)
  dx <- dim(x)
  ddx <- length(dx)
  dn <- if(length(idxs))
          lapply(1:ddx,function(d){
            paste(inms[d],seq.int(dx[d]),sep="=")
          })
        else NULL
  attributes(x) <- NULL
  dim(x) <- dx
  dimnames(x) <- dn
  x
}

print.tensor <- function(x,...){
  
  idxs <- indices(x)
  inms <- index.names(x)
  dx <- dim(x)
  ddx <- length(dx)

  if(length(inms)){
    dn <- lapply(1:ddx,function(d){
      paste(inms[d],seq.int(dx[d]),sep="=")
    })
  }
  else {
    dn <- NULL    
  }
  attributes(x) <- NULL
  dim(x) <- dx
  dimnames(x) <- dn
  
  cat("\nTensor, rank:",ddx,"- dim:",dx)
  if(length(idxs))cat(" - indices:",paste(idxs,collapse=""))
  else cat(" - no indices")
  cat("\n\n")
  print.default(x)
}

print.tensor_skeleton <- function(x,...){
  
  idxs <- indices(x)
  
  cat("\nTensor skeleton, rank:",length(idxs))
  if(length(idxs)) cat(" - indices:",paste(idxs,collapse=""))
  cat("\n<no data>\n")
}

sum.tensor <- function(...,over=NULL){
  
  x <- ..1
  idx <- indices(x)
  if(length(over))
    cl <- index.names(x) %in% over
  else
    cl <- index.types(x)=="+"
  if(!any(cl)) return(x)
  
  x <- collapse.array(unclass(x),which(cl))
  lx <- length(x)
  if(is.array(x))
    structure(x,indices=idx[!cl],class="tensor")
  else if(lx>1)
    structure(x,dim=lx,indices=idx[!cl],class="tensor")
  else x
}

tperm <- function(x,iperm){
  i.names <- index.names(x)
  if(!length(i.names))stop("cannot tperm unnamed tensor")
  perm <- match(iperm,i.names)
  if(any(is.na(perm)))stop("invalid iperm argument")
  
  structure(aperm(x,perm),
            indices=indices(x)[perm],
            class="tensor")
}


fidx <- function(what,where)which(what==where)

perhaps.diag <- function(x){
  
  dd <- seq_along(dim(x))
  ixs <- attr(x,"indices")
  inms <- index.names(x)
  itps <- index.types(x)
  if(!anyDuplicated(inms)) return(x)
  
  dup <- duplicated(inms)
  dup.names <- inms[dup]
  
  uq.dd <- dd[!dup]
  uq.names <- inms[!dup]
  uq.types <- itps[!dup]
  
  has.dup <- sapply(uq.names,`%in%`,dup.names)
  has.dup.dd <- uq.dd[has.dup]
  diag.names <- uq.names[has.dup]
  
  diag.dims <- lapply(diag.names,fidx,inms)
  diag.types <- lapply(diag.dims,elems,itps)
  diag.types <- sapply(diag.types,combine.type)

  rem.names <- uq.names[!has.dup]
  rem.types <- uq.types[!has.dup]
  
  x <- array.diag(x,diag.dims=diag.dims)
  
  res.names <- c(rem.names,diag.names)
  res.types <- c(rem.types,diag.types)
  
  
  
  if(length(x)==1) x
  else {
    
    if(!length(dim(x))) dim(x) <- length(x)
    
    r.ixs <- paste0(res.types,res.names)
    
    perm <- match(uq.names,res.names)
    
    x <- aperm(x,perm)
    r.ixs <- r.ixs[perm]
    
    structure(x,
              indices=r.ixs,
              class="tensor")
  }
}

which.in <- function(what,where)which(what%in%where)
fill_tensor_tensor <- function(x,ix,value){
  
  it <- split.indices(ix)
  inames <- it$names
  itypes <- it$types
  
  if(!anyDuplicated(inames)){
    
    x <- add_indices(x,ix,override=TRUE)
    
    xnames <- index.names(x)
    vnames <- index.names(value)
    if(length(xnames)<length(vnames)) stop("dimensions do not match")
    xvmatch <- match(xnames,vnames)
    
    if(length(xnames)==length(vnames)){
      if(any(is.na(xvmatch))) stop("indices do not match")
      y <- aperm(value,xvmatch)
      if(any(dim(y)!=dim(x))) stop("dimensions do not match")
      dimnames(y) <- NULL
      structure(y,
                indices=indices(x),
                class= "tensor")
    }
    else {
      
      ismatched <- is.finite(xvmatch)
      vperm <- xvmatch[ismatched]
      
      dx <- dim(x)
      dm <- seq_along(dx)
      if(!all(dx[ismatched]==dim(value)[vperm])) stop("dimensions do not match")
      
      xperm <- c(dm[ismatched],dm[!ismatched])
      y <- aperm(x,xperm)
      y[] <- aperm(value,vperm)
      
      inv.xperm <- dm
      inv.xperm[xperm] <- dm
      y <- aperm(y,inv.xperm)
      
      # debugging:
      stopifnot(all(dx==dim(y)))
      dimnames(y) <- NULL
      structure(y,
                indices=attr(x,"indices"),
                class= "tensor")
    }
    
  }
  else {
    
    dup <- duplicated(inames)
    
    dup.names <- inames[dup]
    uq.names <- inames[!dup]
    
    has.dup <- sapply(uq.names,`%in%`,dup.names)
    dup.names <- uq.names[has.dup]
    
    diag.dims <- lapply(dup.names,fidx,inames)
        
    vnames <- index.names(value)
    if(length(uq.names)<length(vnames)) stop("dimensions do not match")
    
    if(length(uq.names)==length(vnames)){
      
      uvmatch <- match(uq.names,vnames)
      if(any(is.na(uvmatch))) stop("indices do not match")
      value <- aperm(value,uvmatch)
    }
    else {
    
      if(!all(vnames %in% uq.names)) stop("indices do not match")
      
      dxu <- dim(x)[!dup]
      is.matched <- uq.names %in% vnames
      value <- array(value,dim=c(dim(value),dxu[!is.matched]))
      vperm <- match(uq.names,c(vnames,uq.names[!is.matched]))
      value <- aperm(value,vperm)
    }
    
    y <- set.array.diag(x,diag.dims,value)
    
    dimnames(y) <- NULL
    structure(y,
              indices=attr(x,"indices"),
              class= "tensor")
  }
  
}


fill_tensor_numeric <- function(x,ix,value){
  
  it <- split.indices(ix)
  inames <- it$names
  itypes <- it$types
  
  if(!anyDuplicated(inames)){
    
    add_indices(array(value,dim=dim(x)),i=ix)
    
  }
  else {
    
    dup <- duplicated(inames)
    
    dup.names <- inames[dup]
    uq.names <- inames[!dup]
    
    has.dup <- sapply(uq.names,`%in%`,dup.names)
    dup.names <- uq.names[has.dup]
    
    diag.dims <- lapply(dup.names,fidx,inames)
    
    y <- set.array.diag(x,diag.dims,value)
    
    dimnames(y) <- NULL
    structure(y,
              indices=attr(x,"indices"),
              class= "tensor")
  }
  
}


fill_tensor_skeleton <- function(x,ix,value){
  
  if(!inherits(value,"tensor")) stop("right hand side must be a tensor")
  
  x <- add_indices(x,ix,override=TRUE)
  
  xnames <- index.names(x)
  vnames <- index.names(value)
  if(length(xnames)!=length(vnames)) stop("dimensions do not match")
  xvmatch <- match(xnames,vnames)
  
  if(any(is.na(xvmatch))) stop("indices do not match")
  
  y <- aperm(value,xvmatch)
  dimnames(y) <- NULL
  structure(y,
            indices=attr(x,"indices"),
            class= "tensor")
}

rowsum.tensor <- function(x, group, reorder=TRUE, ...){
  
    ix <- indices(x)
    x <- rowsum.array(x,group,reorder,...)
    structure(x,indices=ix,class="tensor")
}
  