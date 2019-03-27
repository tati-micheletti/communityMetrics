distanceMatrix<- function(sim){
nSpp<-dim(sim$traits)[1]
nTraits <- dim(sim$traits)[2]-1 #column 1 is species names
d <- matrix(NA, nSpp, nSpp)
for (i in 1:nSpp){
  for (j in 1:nSpp){
    d[i,j] <- sum(sim$traits[i,-1]!=sim$traits[j,-1])
  }
}
d <- d/nTraits 
##estracting values for the lower triangular   
dist.m <- tril(d, -1)
colnames(dist.m)<- sim$traits[,1]
rownames(dist.m)<- sim$traits[,1]
sim$dist.m <- dist.m
return(invisible(sim))
}