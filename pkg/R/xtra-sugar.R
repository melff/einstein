"%[%" <- function(x,i,...)tensor(x)[i,...]

"%[%<-" <- function(x,i,...,value){
  tensor(x)[i,...]<-value
  x
}