#' SIR simulation
#' 
#' Description
#' 
#' @param n Number of individuals.
#' @param beta Infectivity rate.
#' @param gamma Recovery rate.
#' @param status 1xn vector of each individuals state: 0=susceptible, 1=infective, 2=recovery
#' @return status 1xn vector of each individuals state after the epidemic. Final epidemic size is sum(status==2)
#' @author Brock Hermans <brock.hermans@@adelaide.edu.au> 
#' @export
#' @note 9 December 2014
#' @examples
#' # Test that the SIR code gives expected outcomes
#' beta <- seq(1,7,1)
#' par(mfrow=c(3,3))
#' status=rep(0,100); status[1]=1
#' for (b in beta){
#' Ne <- c()
#' for (j in 1:10000){
#' tmp=sir(100,b,1,status)
#' Ne <- c(Ne,sum(tmp==2))
#' }
#' hist(Ne,main=paste("Histogram for Ne for beta=",b))
#' }
sir <- function(n,beta,gamma,status){
  tmp <- seq(1,n)
  time <- c()
  beta <- beta/(n-1)
  
  while (sum(status==1)>0 & sum(status==0)>0){
    num_infec <- sum(status==1)
    num_susc <- sum(status==0)
    infec_rate <- num_susc*beta*num_infec
    recov_rate <- num_infec*gamma
    time <- c(time, sum(time)+rexp(1,num_infec*gamma+num_susc*beta))
    u <- runif(1)
    
    if (u<(infec_rate/(infec_rate+recov_rate))){
      if (sum(status==0)==1){
        infection <- which(status==0)
      }
      else {
        infection <- sample(tmp[status==0],1)
      }
      status[infection] <- status[infection]+1
    }
    else {
      if (sum(status==1)==1){
        recovery <- which(status==1)
      }
      else {
        recovery <- sample(tmp[status==1],1)
      }
      status[recovery] <- status[recovery]+1
    }
    if (sum(status==0)==0){
      status <- rep(2,n)
    }
  }
  if (sum(status==1)>0){
    print("error")
  }
  
  return(status)
}