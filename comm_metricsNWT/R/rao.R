#RAO FUNCTION WORKING BUT SLOW!!!! 
raoEntropy <- function(D, p){
    #p a vector of species proportional abundance, sums to 1
    #D a SxS matrix of distances in trait space
    p2 <- p^2
    S <- dim(D)[1]
    Rao <- 0
     for (i in 1:(S-1))
     for (j in (i+1):S)
     Rao <- Rao + D[i,j]*p2
    return(Rao)
}

raoEntropyM <- function(D, p){
  #p a vector of species proportional abundance, sums to 1
  #D a SxS matrix of distances in trait space
  res <- NA
  if (!is.na(p[1])) {
    res <- p %*% D %*% p
    res <- as.numeric(res)
  }
  return(res)
}

raoEntropyC <- function(D, p){
   rf <- function(x){
     res <- NA
     if (!is.na(x[1])) {
       res <- x %*% D %*% x
       res <- as.numeric(res)
     }
     return(res)
   }
   R <- calc(p, rf)
   return(R)
}

#rao <- function(x, D){
#  x <- x^2
 # S<- dim(D)[1]
 # for (i in 1:(S-1))
  #  for (j in (i+1):S)
   #   rao<-  D[i,j] * x
#}