##
Rcpp::sourceCpp('R/finalSize_mle.cpp')
Rcpp::sourceCpp('R/sim_epidemic.cpp')
source('~/Desktop/DiseaseNetworkInference/diseaseNetworks/R/mleBetaNe.R', echo=TRUE)
source('~/Desktop/DiseaseNetworkInference/diseaseNetworks/R/finalEpidemicSizeDistribution.R', echo=TRUE)

##
status=rep(0,100)
status[1]=1
n = 99
df = as.matrix(expand.grid(0:n,0:n))
df = cbind(df,0)
n=100

##
beta <- seq(1,7,0.1)
betaAold <- c()
timeOld <- 0
betaAnew <- c()
timeNew <- 0
for (b in beta){
  Ne <- c()
  for (i in 1:10){
    tmp <- sim_epidemic(edges=df,node_states=status,gamma=1,Beta=b,trace=0)
    Ne <- c(Ne,sum(tmp==2))
  }
  s1 <- Sys.time()
  betaAold <- c(betaAold,mleBetaNe(Ne))
  f1 <- Sys.time()
  s2 <- Sys.time()
  betaAnew <- c(betaAnew,mle(Ne,beta))
  f2 <- Sys.time()
  timeOld <- timeOld+(f1-s1)
  timeNew <- timeNew+(f2-s2)
}
plot(betaAold,betaAnew,xlab=paste("time=",timeOld),ylab=paste("time=",timeNew))
lines(c(0,8),c(0,8))
