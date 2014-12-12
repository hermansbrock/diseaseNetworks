library(gridExtra)
library(ggplot2)
source('~/Desktop/DiseaseNetworkInference/diseaseNetworks/R/ABC/sir.R', echo=TRUE)
source('~/Desktop/DiseaseNetworkInference/diseaseNetworks/R/ABC/abc_smc.R', echo=TRUE)

n <- 100
status=rep(0,n)
status[sample(n,1)]=1
betaA <- c()

beta <- seq(1,7,1)
N <- 500
T2 <- 4
epsilon <- matrix(c(20,15,11,9,30,25,17,13,100,60,56,36,140,100,80,50,200,161,110,50,210,190,130,50,230,170,120,40),byrow=T,ncol=T2)
#T2 <- 9
#epsilon <- matrix(c(230,200,170,120,90,70,40,23,6),byrow=T,ncol=T2)
out <- list()

for (e in 1:1){
  out[[e]] <- data.frame(matrix(nrow=N))
  for (b in beta){
    D0 <- c()
    for (i in 1:10){
      tmp <- sir(n,b,1,status)
      D0 <- c(D0,sum(tmp==2))
    }
    D0 <- sort(D0)
    smc <- abc_smc(D0,epsilon[b,],N,T2,'Gaussian')
    out[[e]] <- cbind(out[[e]],smc$theta[T2,])
  }
  out[[e]] <- out[[e]][,-1]
  colnames(out[[e]]) <- c("beta1","beta2","beta3","beta4","beta5","beta6","beta7")
  
  
  m1 <- ggplot(out[[e]],aes(x=beta1))
  m2 <- ggplot(out[[e]],aes(x=beta2))
  m3 <- ggplot(out[[e]],aes(x=beta3))
  m4 <- ggplot(out[[e]],aes(x=beta4))
  m5 <- ggplot(out[[e]],aes(x=beta5))
  m6 <- ggplot(out[[e]],aes(x=beta6))
  m7 <- ggplot(out[[e]],aes(x=beta7))
  tmp <- data.frame(matrix(seq(1,7,1),nrow=7))
  tmp2 <- c()
  for (i in 1:7){
    tmp2 <- c(tmp2,mean(out[[e]][,i]))
  }
  tmp <- data.frame(cbind(tmp,tmp2))
  colnames(tmp) <- c("X1","X2")
  m8 <- ggplot(tmp,aes(x=X1,y=X2))
  grid.arrange(
    m1 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(0,8)+ggtitle(paste("beta=",beta[1]))+geom_density(),
    m2 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(0,8)+ggtitle(paste("beta=",beta[2]))+geom_density(),
    m3 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(0,8)+ggtitle(paste("beta=",beta[3]))+geom_density(),
    m4 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(0,8)+ggtitle(paste("beta=",beta[4]))+geom_density(),
    m5 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(0,8)+ggtitle(paste("beta=",beta[5]))+geom_density(),
    m6 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(0,8)+ggtitle(paste("beta=",beta[6]))+geom_density(),
    m7 + geom_histogram(binwidth=1,fill="white",color="black")+xlim(0,8)+ggtitle(paste("beta=",beta[7]))+geom_density(),
    m8 + geom_point()+geom_abline(intercept=0,slope=1)+theme_bw()+xlim(0,8)+ylim(0,8)+xlab("beta")+ylab("beta_ABC"),
    ncol=3,main="ABC SMC"
  )
  
}