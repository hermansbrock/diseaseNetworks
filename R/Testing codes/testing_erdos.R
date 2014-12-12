load("trainSetErdos.RData")
library(MASS)
library(ggplot2)
library(gridExtra)

rownames(trainSetA) <- seq(1,nrow(trainSetA),1)
modOriginal <- lm(beta~betaAux*(meanCon+infecCon),data=trainSetA)

mod <- lm(beta~betaAux*(meanCon+infecCon+pairs+triples+triangles),data=trainSetA)
mod <- stepAIC(mod)
summary(mod)
modNew <- update(mod,.~.-betaAux:infecCon)

summary(modOriginal)
summary(mod)

#Testing for original model
outOriginal <- data.frame(matrix(rep(seq(1,7,0.1),10),ncol=1))
for (i in unique(trainSetA$percent)){
  ind <- which(trainSetA$percent==i)
  betaA <- predict(modOriginal,newdata=trainSetA[ind,])
  outOriginal <- cbind(outOriginal,betaA)
}
rownames(outOriginal) <- seq(1,610,1)
colnames(outOriginal) <- c("beta",seq(5,95,5)/100)

m1 <- ggplot(outOriginal,aes(beta,outOriginal[,3])); m2 <- ggplot(outOriginal,aes(beta,outOriginal[,5])); m3 <- ggplot(outOriginal,aes(beta,outOriginal[,7])); m4 <- ggplot(outOriginal,aes(beta,outOriginal[,9])); m5 <- ggplot(outOriginal,aes(beta,outOriginal[,11])); m6 <- ggplot(outOriginal,aes(beta,outOriginal[,13])); m7 <- ggplot(outOriginal,aes(beta,outOriginal[,15])); m8 <- ggplot(outOriginal,aes(beta,outOriginal[,17])); m9 <- ggplot(outOriginal,aes(beta,outOriginal[,19]))
grid.arrange(
  m1+geom_point()+ggtitle(paste("Plot with percent=",0.1))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m2+geom_point()+ggtitle(paste("Plot with percent=",0.2))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m3+geom_point()+ggtitle(paste("Plot with percent=",0.3))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m4+geom_point()+ggtitle(paste("Plot with percent=",0.4))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m5+geom_point()+ggtitle(paste("Plot with percent=",0.5))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m6+geom_point()+ggtitle(paste("Plot with percent=",0.6))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m7+geom_point()+ggtitle(paste("Plot with percent=",0.7))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m8+geom_point()+ggtitle(paste("Plot with percent=",0.8))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m9+geom_point()+ggtitle(paste("Plot with percent=",0.9))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),  
  ncol=3,main="Original model")



#Testing for new model
outNew <- data.frame(matrix(rep(seq(1,7,0.1),10),ncol=1))
for (i in unique(trainSetA$percent)){
  ind <- which(trainSetA$percent==i)
  betaA <- predict(modNew,newdata=trainSetA[ind,])
  outNew <- cbind(outNew,betaA)
}
rownames(outNew) <- seq(1,610,1)
colnames(outNew) <- c("beta",seq(5,95,5)/100)

m1 <- ggplot(outNew,aes(beta,outNew[,3])); m2 <- ggplot(outNew,aes(beta,outNew[,5])); m3 <- ggplot(outNew,aes(beta,outNew[,7])); m4 <- ggplot(outNew,aes(beta,outNew[,9])); m5 <- ggplot(outNew,aes(beta,outNew[,11])); m6 <- ggplot(outNew,aes(beta,outNew[,13])); m7 <- ggplot(outNew,aes(beta,outNew[,15])); m8 <- ggplot(outNew,aes(beta,outNew[,17])); m9 <- ggplot(outNew,aes(beta,outNew[,19]))
grid.arrange(
  m1+geom_point()+ggtitle(paste("Plot with percent=",0.1))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m2+geom_point()+ggtitle(paste("Plot with percent=",0.2))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m3+geom_point()+ggtitle(paste("Plot with percent=",0.3))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m4+geom_point()+ggtitle(paste("Plot with percent=",0.4))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m5+geom_point()+ggtitle(paste("Plot with percent=",0.5))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m6+geom_point()+ggtitle(paste("Plot with percent=",0.6))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m7+geom_point()+ggtitle(paste("Plot with percent=",0.7))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m8+geom_point()+ggtitle(paste("Plot with percent=",0.8))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),
  m9+geom_point()+ggtitle(paste("Plot with percent=",0.9))+theme_bw()+geom_abline(intercept=0,slope=1)+xlab("beta")+ylab("betaLm")+ylim(0,8),  
  ncol=3,main="New model")