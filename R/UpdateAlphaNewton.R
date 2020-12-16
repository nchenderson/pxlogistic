UpdateAlphaNewton <- function(alpha.old, x.theta, y, n.trials, init=FALSE) {
  if(init){
    alpha.old <- 0
    loop <- TRUE
  }
  else loop <- FALSE
  
  while(TRUE){
     ept <- plogis(alpha.old*x.theta)
     dept <- ept*(1-ept)
     
     f.st <- sum(x.theta*(y - n.trials*ept))
     if(loop==FALSE | abs(f.st) < 1e-7) break
     
     f.nd <- sum(n.trials*x.theta*x.theta*dept)
     alpha.new <- as.numeric(alpha.old + f.st/f.nd)
     alpha.old <- alpha.new
  }
  ans <- alpha.new
  return(ans)
}
