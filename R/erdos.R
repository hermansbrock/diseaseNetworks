erdos <- function(n,p){
  nodes <- matrix(runif(n*2),ncol=2)
  edges <- matrix(nrow=1,ncol=2)
  A <- matrix(0,ncol=n,nrow=n)
  for (i in 1:n){
    for (j in 1:n){
      u <- runif(1)
      if (u<p){
        edges <- rbind(edges,c(i,j))
        edges <- rbind(edges,c(j,i))
        A[i,j] <- 1
        A[j,i] <- 1
      }
    }
  }
  x <- list()
  x$nodes <- nodes
  x$edges <- edges
  x$A <- A
  x$percent <- (nrow(edges)/2)/9900
  return(x)
}