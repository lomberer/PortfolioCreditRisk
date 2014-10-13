
N <- 100 #number of Loans
M <- 1000 #number of scenario
#Z <- -1.5 #Credit Factor

PosNames <- paste("Loan",1:N,sep="")
rho <- rep(0.3,nrow=N) #correlation with Z
LGD <- rep(0.5,nrow=N) #% of Loss
PD <- rep(0.01,nrow=N) #Probability of Default
EAD <- rep(100,nrow=N)  #Exposure at Default

port <- data.frame(PosNames,PD,LGD,EAD,rho,stringsAsFactors = F)
head(port)

simMC <- function(port,M){
  N <- nrow(port)
  Z <- (rnorm(M,mean=0))
  LossMC <- matrix(,M,1)
  for (m in 1:M){
    e <- rnorm(N,mean=0,sd=1)
    V <- rho*Z[m] + sqrt(1-rho^2)*e
    #default_treshold <- (qnorm(PD)-Z*t(rho))/sqrt(1-rho^2)
    
    default_flag <- V < qnorm(PD)
    LossMC[m] <- sum(default_flag * port$EAD * port$LGD)
  }
  return (cbind(LossMC,1/M))
}

simIS <- function(port,M,mu){
  N <- nrow(port)
  if (exists("SEED")) {set.seed(SEED)
  } else set.seed(12345)

  Z <- (rnorm(M,mean=mu))
  LossMC <- matrix(,M,1)
  for (m in 1:M){
    default_treshold <- (qnorm(PD)-Z[m]*rho)/sqrt(1-rho^2)
    e_matrix <- matrix(rnorm(N,mean=0,sd=1),1,N)
    default_flag <- e_matrix <default_treshold
    LossMC[m] <- sum(default_flag * port$EAD * port$LGD)
  }
  ISweight <- exp(-mu*Z+mu^2/2)/M
  return (cbind(LossMC,ISweight))
}

simRes <- function(simLoss,qt=c(0.8,0.9,0.95,0.99,0.999,0.9999)){
  return(wtd.quantile(simLoss[,1],weights=simLoss[,2],probs=qt))
}

wtable <- function(simLoss){
  return(wtd.table(simLoss[,1],weights=simLoss[,2]))
}

naiveMC <- function(port,M) {
  lossMC <- matrix(,M)
  N <- nrow(port)
  for (m in 1:M){
    loss_m=0
    for (n in 1:N){
      Z <- rnorm(1)
      e <- rnorm(1)
      w <- port$rho[n]
      if (w*Z+sqrt(1-w^2)*e < qnorm(port$PD[n])){
        loss_m = loss_m + port$LGD[n] *port$EAD[n]
      }
    }
    lossMC[m]=loss_m
    }
    return (cbind(lossMC,1/M))
  }

simMC <- function(port,M,mu){
  N <- nrow(port)
  Z <- (rnorm(M,mean=mu))
  LossMC <- matrix(,M,1)
  for (m in 1:M){
    e <- rnorm(N,mean=0,sd=1)
    V <- rho*Z[m] + sqrt(1-rho^2)*e
    #default_treshold <- (qnorm(PD)-Z*t(rho))/sqrt(1-rho^2)
    default_flag <- V < qnorm(PD)
    LossMC[m] <- sum(default_flag * port$EAD * port$LGD)
  }
  weights <- exp(-mu*Z+mu^2/2)/M
  return (cbind(LossMC,weights))
}





"
M = 9999
LossMC <- simMC(port,M)
LossIS <- simIS(port,M,-1.5)
simRes(LossMC)
simRes(LossIS)
"