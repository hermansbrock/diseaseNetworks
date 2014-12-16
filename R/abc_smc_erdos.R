abc_smc_erdos <- function(D0,epsilon,N,T2,METHOD){
  A_percent <- matrix(nrow=T2,ncol=N)
  W_percent <- matrix(nrow=T2,ncol=N)
  
  A_beta <- matrix(nrow=T2,ncol=N)
  W_beta <- matrix(nrow=T2,ncol=N)
  
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