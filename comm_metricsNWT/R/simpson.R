simpson <- function(p){
  D <-p^2 # Simpson dominance index per species
  D_in <- 1- sum(D)  # Simpson Diversity index
  
  return(D_in)
}  