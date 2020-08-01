UpdateAlphaNewton <- function(alpha.old, x.theta, y, init=FALSE){
  if(init){
    alpha.old <- 0
    loop <- TRUE
  }
  else loop <- FALSE
  
  while(TRUE){
     ept <- expit(alpha.old*x.theta)
     dept <- ept*(1-ept)
     f.st <- sum(x.theta*(y - n.trials*ept))
     f.nd <- sum(n.trials*x.theta*x.theta*dept)
     alpha.new <- as.numeric(alpha.old + f.st/f.nd)
     if(loop==FALSE | norm(alpha.new-alpha.old, "2") < 1e-6) break
     alpha.old <- alpha.new
  }
  ans <- alpha.new
  return(ans)
}

expit <- function(x) {
   1/(1 + exp(-x))
}