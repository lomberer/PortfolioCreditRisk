d = 10 # 1...j...d
m = 1000 # 1...k...m
x = 1000 #loss treshold
p = rep(NA,m)
c = rep(NA,m)
for (k in 1:m){
  p[k] = 0.01*(1+sin(16*pi*k/m)) #PD
  c[k] = ceiling(5*k/m)^2
}
a = matrix(runif(d*m,min=0,max=1/sqrt(d)),nrow=m) #1000 × 10
b = sqrt(1-rowSums(a^2)) #1000



getp_theta <- function(theta,p,c){
  return(p*exp(theta*c)/(1+p*(exp(theta*c)-1)))
}

getpsi <- function(theta,p,c){
  psi=0
  for (k in 1:m){
    psi = psi + log(1+p[k]*(exp(theta*c[k])-1))
  }
  return(psi)
}


getlr <- function(theta,L,p,c){
  return(exp(-theta*L+getpsi(theta,p,c)))
}

#\psi.prime
getdpsi <- function(theta,p,c){
  dpsi = 0
  for (k in 1:length(p)){
    dpsi = dpsi+ p[k]*c[k]*exp(theta*c[k])/(1+p[k]*(exp(theta*c[k])-1))
  }
  return(dpsi)
}

gettheta_x <- function(x,p_Z,c){
  oz <- function(theta){
    dpsi = 0 
    for (k in 1:length(p_Z)){
      dpsi = dpsi+ p_Z[k]*c[k]*exp(theta*c[k])/(1+p_Z[k]*(exp(theta*c[k])-1))
    }
    return(dpsi-x)
  }
  if (sum(p_Z*c)>x){theta_x = 0}
  else {
    theta_x <- uniroot(oz,c(0,5))$root
  }
  return(theta_x)
}

getJ_x <- function(z){
  #Z is vector of d
  p_Z <- pnorm((a%*%Z+qnorm(p))/b)
  oz <- function(theta){
    dpsi = 0 
    for (k in 1:length(p_Z)){
      dpsi = dpsi+ p_Z[k]*c[k]*exp(theta*c[k])/(1+p_Z[k]*(exp(theta*c[k])-1))
    }
    return(dpsi-x)
  }
  if (sum(p_Z*c)>x){theta_x = 0}
  else {
    theta_x <- uniroot(oz,c(0,2))$root
  }
  psi <- getpsi(theta_x,p_Z,c)
  return(-1*(-theta_x*x+getpsi(theta_x,p_Z,c)-t(Z)%*%Z))
}
#Calibrate \mu by maximizing J(x)
mu <- optim(runif(d, 0, 1),getJ_x,method = "L-BFGS-B",lower=rep(0,d),upper=rep(2,d))$par

mu=rep(0,d)
#try to reduct variance only by the first part.
Glassman_MC <- function(){
  #Step 1
  Z = rnorm(d,mean=mu,sd=1)
  #Step 2
  p_Z <- pnorm((a%*%Z+qnorm(p))/b)
  theta_x <- gettheta_x(x,p_Z,c)
  #print(theta_x)
  p_theta <- getp_theta(theta_x,p_Z,c)
  #Step 3
  X <- a%*%Z+b*rnorm(m)
  Y <- rnorm(m)>qnorm(1-p_theta)
  L <- sum(Y*c)
  #print(L)
  weight=exp(-theta_x*L+getpsi(theta_x,p_Z,c))#weight is too low. 
  res <- c(L,weight)
  names(res) <- c("LossAmount","LossProb")
  return(res)
}

N = 100
resv=data.frame(Loss=rep(NA,N),Weights=rep(NA,N))
for (i in 1:N){
  resv[i,1] = Glassman_MC()[1]
  resv[i,2] = Glassman_MC()[2]
}
lossquantile(resv)
