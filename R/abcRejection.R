#' ABC Rejection for SIR
#' 
#' This function implements SIR Rejection for the SIR model based on final epidemic size.
#' 
#' @param D0 Original dataset.
#' @param epsilon tolerance level.
#' @param N Number of samples from the posterior that we need to sample.
#' @return theta N samples from the posterior distribution.
#' @author Brock Hermans <brock.hermans@@adelaide.edu.au
#' @export
#' @note 1 December 2014
#' @examples
#' library(gridExtra); library(ggplot2)
#' beta <- seq(1,7,1); n <- 100; epsilon <- 20; status=rep(0,n); status[sample(n,1)]=1
#' par(mfrow=c(3,3))
#' out <- list()
#' N <- 500
#' e=epsilon
#' out[[1]] <- data.frame(matrix(nrow=N))
#' betaA <- c()
#' for (b in beta){
#'  D0 <- c()
#'  for (i in 1:10){
#'   tmp <- sir(n,b,1,status)
#'   D0 <- c(D0,sum(tmp==2))}
#'  D0 <- sort(D0)
#'  theta <- abcRejection10(D0,e,N)
#'  out[[1]] <- cbind(out[[1]],theta)
#'  betaA <- c(betaA,mean(theta))}
#' out[[1]] <- out[[1]][,-1]
#' colnames(out[[1]]) <- c("beta1","beta2","beta3","beta4","beta5","beta6","beta7")
#' m1 <- ggplot(out[[1]],aes(x=beta1)); m2 <- ggplot(out[[1]],aes(x=beta2)); m3 <- ggplot(out[[1]],aes(x=beta3));  m4 <- ggplot(out[[1]],aes(x=beta4)); m5 <- ggplot(out[[1]],aes(x=beta5)); m6 <- ggplot(out[[1]],aes(x=beta6)); m7 <- ggplot(out[[1]],aes(x=beta7))
#' tmp <- data.frame(matrix(seq(1,7,1),nrow=7))
#' tmp2 <- c()
#' for (i in 1:7){
#'  tmp2 <- c(tmp2,mean(out[[1]][,i]))}
#' tmp <- data.frame(cbind(tmp,tmp2))
#' colnames(tmp) <- c("X1","X2")
#' m8 <- ggplot(tmp,aes(x=X1,y=X2))
#' grid.arrange(
#' m1 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(1,7)+ggtitle(paste("beta=",beta[1]))+geom_density(),
#' m2 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(1,7)+ggtitle(paste("beta=",beta[2]))+geom_density(),
#' m3 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(1,7)+ggtitle(paste("beta=",beta[3]))+geom_density(),
#' m4 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(1,7)+ggtitle(paste("beta=",beta[4]))+geom_density(),
#' m5 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(1,7)+ggtitle(paste("beta=",beta[5]))+geom_density(),
#' m6 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(1,7)+ggtitle(paste("beta=",beta[6]))+geom_density(),
#' m7 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(1,7)+ggtitle(paste("beta=",beta[7]))+geom_density(),
#' m8 + geom_point()+geom_abline(intercept=0,slope=1)+theme_bw()+xlab("beta")+ylab("beta_ABC"),
#' ncol=3,main=paste("epsilon=",e))
abcRejection10 <- function(D0,epsilon,N){
  n=100
  status=rep(0,n)
  status[1]=1
  
  theta <- c()
  
  for (i in 1:N){
    print(i)
    ind=F
    while (ind==F){
      thetaStar <- runif(1,1,7)
      NeStar <- c()
      for (k in 1:10){
        tmp = sir(n,thetaStar,1,status)
        NeStar <- c(NeStar,sum(tmp==2))
      }
      NeStar <- sort(NeStar)
      if (dist(rbind(NeStar,D0))<=epsilon){
        theta <- c(theta,thetaStar) 
        ind=T
      }
    }
  }
  return(theta)
}
