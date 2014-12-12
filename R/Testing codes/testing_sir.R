beta <- seq(1,7,1)
par(mfrow=c(3,3))
status=rep(0,100); status[1]=1
for (b in beta){
  Ne <- c()
  for (j in 1:10000){
    tmp=sir(100,b,1,status)
    Ne <- c(Ne,sum(tmp==2))
  }
    hist(Ne,main=paste("Histogram for Ne for beta=",b))
}

