m = 1000 #number of loans
p = rep(NA,m) #PD
c = rep(NA,m) #EAD*LGD
for (k in 1:m){
  p[k] = 0.01*(1+sin(16*pi*k/m)) #PD
  c[k] = ceiling(5*k/m)^2
}
p = rep(0.01,m)

d = 10 #10-factor model
a = matrix(rep(0.1,m*d),nrow=d)
#a = matrix(runif(m*d,min=0,max=1/sqrt(d)),nrow=d)#matrix of 10x1000
b = sqrt(1-colSums(a^2)) #vector of 1000x1
x = 1000 #loss treshold

#==================
#eq(2)
getp_Z <- function(p,Z,a,b){
  return(pnorm((colSums(a*Z)+qnorm(p))/b))
  "  m <- length(p)
  p_Z <- rep(NA,m)
  for (k in 1:m){
  p_Z[k]=pnorm((sum(a[,k]*Z[,k])+qnorm(p[k]))/b[k])
  }
  return(p_Z)
  "}



getpsi <- function(theta, p,c){
  psi = 0
  for (k in 1:length(p)){
    psi = psi + log(1+p[k]*(exp(theta*c[k])-1))
  }
  return(psi)
}

#likelyhood ratio in eq(5)
getlr <- function(theta, p,c,L){
  return(exp(-theta*L+getpsi(theta,p,c)))
}
#\psi.prime
getdpsi <- function(theta,p_Z,c){
  dpsi = 0
  for (k in 1:length(p_Z)){
    dpsi = dpsi+ p_Z[k]*c[k]*exp(theta*c[k])/(1+p_Z[k]*(exp(theta*c[k])-1))
  }
  return(dpsi)
}
#solve for shifted mean = x, return theta_x

gettheta_x <- function(x,p_Z,c){
  oz <- function(theta){
    return(getdpsi(theta,p_Z,c)-x)
  }
  theta_x = uniroot(oz,c(0,3))$root
  return(theta_x)
}
#$p_{k,\theta} in eq(4)
getp_theta <- function(theta,p,c){
  return(p*exp(theta_x*c)/(1+p*(exp(theta_x*c)-1)))
}
#pt <- p_theta(tx,p,c)


#=====================
#Step 1
Z <- matrix(rnorm(d*m),ncol=m)
# Y = colSums(a*Z)+b*rnorm(m,0,1) < qnorm(p)#condition of default
#Step 2
p_Z <- getp_Z(p,Z,a,b)
theta_x <- gettheta_x(x,p_Z,c)
#theta_x <- theta_x(x,p,c)

#Step 3
p_theta <- getp_theta(theta_x,p_Z,c)
Y = colSums(a*Z)+b*rnorm(m,0,1) < qnorm(p_theta)

#Step 4
L = sum(c*Y)
res = (L>x)*exp(-theta_x*L+getpsi(theta_x,p_theta,c))



N=1000
LossMC=matrix(rep(c(1,2),N),nrow=N,byrow=T)
for (i in 1:N){
  LossMC[i,]=GMC()
}

#——----------------------
F_x <- function(theta_x,x,p,c){
  return(-theta_x*x + getpsi(theta_x,p,c))
}



