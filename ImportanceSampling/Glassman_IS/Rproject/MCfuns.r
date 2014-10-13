PortLossMC <- function(port,N,M,seed=12345){
  set.seed(12345)
  e_matrix <- matrix(rnorm(N*M,mean=0,sd=1),M,N)
  default_treshold <- (qnorm(PD)-rho*Z)/sqrt(1-rho^2)
  default_flag <- e_matrix <default_treshold
  
  LossMC <- default_flag * port$EAD * port$LGD
  PortLoss <- rowSums(LossMC)
  return (PortLoss)
}
