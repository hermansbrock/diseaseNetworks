#' ABC SMC for SIR model
#' 
#' Implements ABC SMC for the SIR model based on final epidemic size.
#' 
#' @param D0 Original dataset.
#' @param epsilon Vector of T2 tolerance levels.
#' @param N Number of points to retain at each population.
#' @param T2 Number of populaitons to generate.
#' @param Method String that determines perturbation kernel. Either Method='Gaussian' or 'Uniform'.  
#' @return A T2xN matrix of sampled points.
#' @return W T2xN matrix of weights
#' @author Brock Hermans <brock.hermans@@adelaide.edu.au>
#' @export
#' @note 8 Decemeber 2014
#' @examples
#' library(gridExtra)
#' library(ggplot2)
#' n <- 100; status=rep(0,n); status[sample(n,1)]=1
#' betaA <- c()
#' beta <- seq(1,7,1); N <- 500; T2 <- 9; epsilon <- matrix(c(230,200,170,120,90,70,40,23,6),byrow=T,ncol=T2)
#' ##
#' out <- list()
#' out[[1]] <- data.frame(matrix(nrow=N))
#' for (b in beta){
#'  D0 <- c()
#'  for (i in 1:10){
#'    tmp <- sir(n,b,1,status)
#'    D0 <- c(D0,sum(tmp==2))}
#'  D0 <- sort(D0)
#'  smc <- abc_smc(D0,epsilon[b,],N,T2,'Gaussian')
#'  out[[1]] <- cbind(out[[1]],smc$theta[T2,])}
#' ##
#' out[[1]] <- out[[1]][,-1]
#' colnames(out[[1]]) <- c("beta1","beta2","beta3","beta4","beta5","beta6","beta7")
#' m1 <- ggplot(out[[e]],aes(x=beta1)); m2 <- ggplot(out[[e]],aes(x=beta2)); m3 <- ggplot(out[[e]],aes(x=beta3)); m4 <- ggplot(out[[e]],aes(x=beta4)); m5 <- ggplot(out[[e]],aes(x=beta5)); m6 <- ggplot(out[[e]],aes(x=beta6)); m7 <- ggplot(out[[e]],aes(x=beta7))
#' tmp <- data.frame(matrix(seq(1,7,1),nrow=7))
#' tmp2 <- c()
#' for (i in 1:7){
#'  tmp2 <- c(tmp2,mean(out[[e]][,i]))}
#'  tmp <- data.frame(cbind(tmp,tmp2))
#'  colnames(tmp) <- c("X1","X2")
#'  m8 <- ggplot(tmp,aes(x=X1,y=X2))
#'  grid.arrange(
#'  m1 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(0,8)+ggtitle(paste("beta=",beta[1]))+geom_density(),
#'  m2 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(0,8)+ggtitle(paste("beta=",beta[2]))+geom_density(),
#'  m3 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(0,8)+ggtitle(paste("beta=",beta[3]))+geom_density(),
#'  m4 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(0,8)+ggtitle(paste("beta=",beta[4]))+geom_density(),
#'  m5 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(0,8)+ggtitle(paste("beta=",beta[5]))+geom_density(),
#'  m6 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(0,8)+ggtitle(paste("beta=",beta[6]))+geom_density(),    
#'  m7 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(0,8)+ggtitle(paste("beta=",beta[7]))+geom_density(),
#'  m8 + geom_point()+geom_abline(intercept=0,slope=1)+theme_bw()+xlim(0,8)+ylim(0,8)+xlab("beta")+ylab("beta_ABC"),
#'  ncol=3,main="ABC SMC")
abc_smc <- function(D0,epsilon,N,T2,Method){
  A <- matrix(nrow=T2,ncol=N)
  W <- matrix(nrow=T2,ncol=N)
  
  n=100
  status=rep(0,n)
  status[1]=1
  
  for (t in 1:T2){
    print(t)
    for (i in 1:N){
      print(i)
      indi <- T
      while (indi){
        thetaStSt <- sim_theta(t,A,W,Method)
          #Prior is never zero as we use a uniform prior
          NeStar <- c()
          for (j in 1:10){
            tmp <- sir(100,thetaStSt,1,status)
            NeStar <- c(NeStar,sum(tmp==2))
          }
          NeStar <- sort(NeStar)
          if (dist(rbind(D0,NeStar))<=epsilon[t]){
            indi <- F
            A[t,i] <- thetaStSt
            W[t,i] <- calculate_weights(thetaStSt,t,A,W,Method)
          }
      }
    }
  }
  out <- list()
  out$theta <- A
  out$weights <- W
  return(out)
} 

pert_kern <- function(t,theta=0,thetaSt,method){
  library(tmvtnorm)
  if (method=='Gaussian'){
    if (theta==0){
      out <- rtmvnorm(1,mean=thetaSt,sigma=1,lower=1,upper=7)
    }
    else {
      out <- dtmvnorm(theta,mean=c(thetaSt),sigma=1,lower=1,upper=7) 
    }
  }
  else if (method=='Uniform'){
    if (theta==0){
      out <- runif(1,max(1,thetaSt-1),min(7,thetaSt+1))
    }
    else {
      out <- dunif(theta,max(1,thetaSt-1),min(7,thetaSt+1))
    }
  }
  else{stop("Invalid method for perturbation")}
  return(out)
}

sim_theta <- function(t,A,W,Method){
  if (t==1){
    out <- runif(1,1,7)
  }
  else{
    thetaSt <- sample(x=A[t-1,],size=1,prob=W[t-1,])
    out <- pert_kern(t,theta=0,thetaSt=thetaSt,method=Method)
  }
  return(out)
}

calculate_weights <- function(thetaStSt,t,A,W,Method){
  if (t==1){
    out <- 1
  }
  else {
    out <- 0
    for (j in 1:ncol(A)){
      out <- out+W[t-1,j]*pert_kern(t,theta=A[t-1,j],thetaSt=thetaStSt,method=Method)
    }
    out <- dunif(thetaStSt,min=1,max=7)/out
  }
  return(out)
}
