library(gridExtra)
library(ggplot2)
source('~/Desktop/DiseaseNetworkInference/diseaseNetworks/R/ABC/abcRejection.R', echo=TRUE)
source('~/Desktop/DiseaseNetworkInference/diseaseNetworks/R/ABC/sir.R', echo=TRUE)

beta <- seq(1,7,1)
n <- 100
epsilon <- 20
status=rep(0,n)
status[sample(n,1)]=1
betaA <- c()

par(mfrow=c(3,3))
out <- list()

N <- 500

for (e in epsilon){
  print(e)
  ind <- which(epsilon==e)
  out[[ind]] <- data.frame(matrix(nrow=N))
  betaA <- c()
  for (b in beta){
    print(b)
    D0 <- c()
    for (i in 1:10){
    tmp <- sir(n,b,1,status)
    D0 <- c(D0,sum(tmp==2))
    }
    D0 <- sort(D0)
    theta <- abcRejection10(D0,e,N)
    out[[ind]] <- cbind(out[[ind]],theta)
    betaA <- c(betaA,mean(theta))
#    hist(theta,main=paste("beta=",b))
  }
  out[[ind]] <- out[[ind]][,-1]
  colnames(out[[ind]]) <- c("beta1","beta2","beta3","beta4","beta5","beta6","beta7")

m1 <- ggplot(out[[ind]],aes(x=beta1))
m2 <- ggplot(out[[ind]],aes(x=beta2))
m3 <- ggplot(out[[ind]],aes(x=beta3))
m4 <- ggplot(out[[ind]],aes(x=beta4))
m5 <- ggplot(out[[ind]],aes(x=beta5))
m6 <- ggplot(out[[ind]],aes(x=beta6))
m7 <- ggplot(out[[ind]],aes(x=beta7))
tmp <- data.frame(matrix(seq(1,7,1),nrow=7))
tmp2 <- c()
for (i in 1:7){
  tmp2 <- c(tmp2,mean(out[[ind]][,i]))
}
tmp <- data.frame(cbind(tmp,tmp2))
colnames(tmp) <- c("X1","X2")
m8 <- ggplot(tmp,aes(x=X1,y=X2))
grid.arrange(
  m1 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(1,7)+ggtitle(paste("beta=",beta[1]))+geom_density(),
  m2 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(1,7)+ggtitle(paste("beta=",beta[2]))+geom_density(),
  m3 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(1,7)+ggtitle(paste("beta=",beta[3]))+geom_density(),
  m4 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(1,7)+ggtitle(paste("beta=",beta[4]))+geom_density(),
  m5 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(1,7)+ggtitle(paste("beta=",beta[5]))+geom_density(),
  m6 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(1,7)+ggtitle(paste("beta=",beta[6]))+geom_density(),
  m7 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(1,7)+ggtitle(paste("beta=",beta[7]))+geom_density(),
  m8 + geom_point()+geom_abline(intercept=0,slope=1)+theme_bw()+xlab("beta")+ylab("beta_ABC"),
  ncol=3,main=paste("epsilon=",e)
  )
}

