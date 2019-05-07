richness <- function(lambdas, scale=6.25){
  #expecting raster stack of predicted mean abundances per ha
  prob0 <- exp(-scale * lambdas)
  p.occ <- 1 - prob0
  return(sum(p.occ)) #we expect sum() to yield a raster  
}