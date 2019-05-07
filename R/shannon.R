
shannon <- function(p){
  #we expect p is a raster stack
    tmp <- p * log(p)
    H <- 1 - sum(tmp)
    return(H)
  }  


#  area<- prod(res(bird.abun))/10000
 # return (area)
#}

