N <- 100 #number of loans.
PosNames <- paste("Loan",1:N,sep="")
rho <- rep(0.3,N) #correlation with macroeconomics
LGD <- rep(0.5,N) 
PD <- rep(0.01,N)
EAD <- rep(100,N)

port <- data.frame(PosNames,PD,LGD,EAD,rho,stringsAsFactors = F)
head(port)

simMC <- function(port,M){
  N <- nrow(port)
  Z <- (rnorm(M,mean=0))
  LossMC <- matrix(,M,1)
  for (m in 1:M){
    e <- rnorm(N,mean=0,sd=1)
    V <- rho*Z[m] + sqrt(1-rho^2)*e
    default_flag <- V < qnorm(PD)
    LossMC[m] <- sum(default_flag * port$EAD * port$LGD)
  }
  return (data.frame(Loss=LossMC,Weights=1/M))
  #return (cbind(LossMC,1/M))
}

simIS <- function(port,M,mu){
  N <- nrow(port)
  Z <- (rnorm(M,mean=mu))
  LossIS <- matrix(,M,1)
  for (m in 1:M){
    e <- rnorm(N,mean=0,sd=1)
    V <- rho*Z[m] + sqrt(1-rho^2)*e
    #default_treshold <- (qnorm(PD)-Z*t(rho))/sqrt(1-rho^2)
    default_flag <- V < qnorm(PD)
    LossIS[m] <- sum(default_flag * port$EAD * port$LGD)
  }
  weights <- exp(-mu*Z+mu^2/2)/M
  return (data.frame(Loss=LossIS,Weights=weights))
}


LossMC <- simMC(port,1000)
LossIS <- simIS(port,1000,-2)
wtd.quantile(LossMC$Loss,weights=LossMC$Weights,probs=c(0.9,0.95,0.99,0.999,0.9999),type='i/n')
wtd.quantile(LossIS$Loss,weights=LossIS$Weights,probs=c(0.9,0.95,0.99,0.999,0.9999),type='i/n')


quantile(LossMC$Loss,c(0.9,0.95,0.99,0.999,0.9999))

library(Hmisc)
library(hydroPSO)
wquantile(LossMC[,1],weights=LossMC[,2],probs=c(0.9,0.95,0.99,0.999),byrow=T)

psi_GMC <- function(theta,p,c){
  return (sum(1+p*(exp(theta * c)-1)))
}