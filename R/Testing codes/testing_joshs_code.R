net <- networkSimulate(1,100)
beta <- seq(1,7,0.1)
status=rep(0,100)
status[1]=1
par(mfrow=c(3,1))
n = 99
df = as.matrix(expand.grid(0:n,0:n))
df = cbind(df,0)
n=100

## Based on 1 outbreak
betaJ1 <- c()
for (i in beta){
  tmp <- sir(100,i,1,status)
  betaJ1 <- c(betaJ1,mleBetaNe(sum(tmp==2)))
}
plot(beta,betaJ1,main="Based on 1 outbreak, sir code",xlim=c(0,8),ylim=c(0,8),xlab=expression(beta),ylab="beta_J")

## Based on 1 outbreak - Brock's code
betaJB <- c()
for (i in beta){
  tmp <- diseaseSimulation(net,i)
  betaJB <- c(betaJB,tmp$Ne)
}
plot(beta,betaJB,main="Based on 1 outbreak, diseaseSimulation",xlim=c(0,8),ylim=c(0,8),xlab=expression(beta),ylab="beta_J")

## Based on 1 outbreak - Jono's code
betaJT <- c()
for (i in beta){
  tmp <- sim_epidemic(edges=df,node_states=status,gamma=1,Beta=i,trace=0)
  betaJT <- c(betaJT,mleBetaNe(sum(tmp==2)))
}
plot(beta,betaJT,main="Based on 1 outbreak, sim_epidemic",xlim=c(0,8),ylim=c(0,8),xlab=expression(beta),ylab="beta_J")


## Based on 10 outbreaks
betaJ10 <- c()
for (i in beta){
  NeTp <- c()
  for (j in 1:10){
    tmp <- sir(100,i,1,status)
    NeTp <- c(NeTp,sum(tmp==2))
  }
  betaJ10 <- c(betaJ10,mleBetaNe(NeTp))
}
plot(beta,betaJ10,main="Based on 10 outbreaks, sir code",xlim=c(0,8),ylim=c(0,8),xlab=expression(beta),ylab="beta_J")

## Based on 10 outbreaks - Brock's code
betaJB10 <- c()
for (i in beta){
  NeTp <- c()
  for (j in 1:10){
    tmp <- diseaseSimulation(net,i)
    NeTp <- c(NeTp,tmp$Ne)
  }
  betaJB10 <- c(betaJB10,mleBetaNe(NeTp))
}
plot(beta,betaJB10,main="Based on 10 outbreak, diseaseSimulation",xlim=c(0,8),ylim=c(0,8),xlab=expression(beta),ylab="beta_J")

## Based on 10 outbreaks - Jono's code
betaJT10 <- c()
for (i in beta){
  NeTp <- c()
  for (j in 1:10){
    tmp <- sim_epidemic(edges=df,node_states=status,gamma=1,Beta=i,trace=0)
    NeTp <- c(NeTp,sum(tmp==2))
  }
  betaJT10 <- c(betaJT10,mleBetaNe(NeTp))
}
plot(beta,betaJT10,main="Based on 10 outbreak, sim_epidemic",xlim=c(0,8),ylim=c(0,8),xlab=expression(beta),ylab="beta_J")

## Based on 1 outbreak repeated 10 times
betaJ10R <- c()
for (i in beta){
  tmp <- sir(100,i,1,status)
  betaJ10R <- c(betaJ10R,mleBetaNe(rep(sum(tmp==2),10)))
}
plot(beta,betaJ10R,main="1 outbreak repeated 10 times, sir code",xlim=c(0,8),ylim=c(0,8),xlab=expression(beta),ylab="beta_J")

## Based on 1 outbreak repeated 10 times - Brock's code
betaJB10R <- c()
for (i in beta){
  tmp <- diseaseSimulation(net,i)
  betaJB10R <- c(betaJB10R,mleBetaNe(rep(tmp$Ne,10)))
}
plot(beta,betaJB10R,main="1 outbreak repeated 10 times, diseaseSimulation",xlim=c(0,8),ylim=c(0,8),xlab=expression(beta),ylab="beta_J")

##Based on 1 outbreak repeated 10 times - Jono's code
betaJT10R <- c()
for (i in beta){
    tmp <- sim_epidemic(edges=df,node_states=status,gamma=1,Beta=i,trace=0)
  betaJT10R <- c(betaJT10R,mleBetaNe(rep(sum(tmp==2),10)))
}
plot(beta,betaJT10R,main="1 outbreak repeated 10 times, sim_epidemic",xlim=c(0,8),ylim=c(0,8),xlab=expression(beta),ylab="beta_J")

