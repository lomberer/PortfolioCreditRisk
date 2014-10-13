N <- 1000 #number of loans.
PosNames <- paste("Loan",1:N,sep="")
beta <- rep(0.3,N) #factor loading
LGD <- rep(0.5,N) 
PD <- rep(0.01,N)
EAD <- rep(100,N)

port <- data.frame(PosNames,PD,LGD,EAD,beta,stringsAsFactors = F)
head(port)

###
simMC <- function(port,M){
  N <- nrow(port)
  Z <- (rnorm(M,mean=0))
  LossMC <- matrix(,M,1)
  for (m in 1:M){
    e <- rnorm(N,mean=0,sd=1)
    V <- beta*Z[m] + sqrt(1-beta^2)*e
    default_flag <- V < qnorm(PD)
    LossMC[m] <- sum(default_flag * port$EAD * port$LGD)
  }
  return (data.frame(Loss=LossMC,Weights=1/M))
}

simIS <- function(port,M,mu){
  N <- nrow(port)
  Z <- (rnorm(M,mean=mu))
  LossIS <- matrix(,M,1)
  for (m in 1:M){
    e <- rnorm(N,mean=0,sd=1)
    V <- beta*Z[m] + sqrt(1-beta^2)*e
    #default_treshold <- (qnorm(PD)-Z*t(beta))/sqrt(1-beta^2)
    default_flag <- V < qnorm(PD)
    LossIS[m] <- sum(default_flag * port$EAD * port$LGD)
  }
  weights <- exp(-mu*Z+mu^2/2)/M
  return (data.frame(Loss=LossIS,Weights=weights))
}

library(Hmisc)
q = c(0.5,0.9,0.95,0.99,0.999,0.9999)
lossplot <- function(LossRes){
  LossAmt <- LossRes$Loss
  LossWgt <- LossRes$Weights
  Ecdf(LossAmt,weights = LossWgt,datadensity='density',xlim=c(0,10000),ylim=c(0.9,1),q=q)
}

lossquantile <- function(LossRes){
  LossAmt <- LossRes$Loss
  LossWgt <- LossRes$Weights
  return(wtd.quantile(LossAmt,weights=LossWgt,probs=q,normwt=T,type='i/n'))
}

LossBench <- simMC(port,100000)
LossMC <- simMC(port,10000)
LossIS <- simIS(port,10000,-3)

Benchmark <- lossquantile(LossBench)
NaiveMC <- lossquantile(LossMC)
ImportanceSampling <- lossquantile(LossIS)
df <- data.frame(Benchmark,NaiveMC,ImportanceSampling)

t(df)

lossplot(LossBench)
lossplot(LossMC)
lossplot(LossIS)