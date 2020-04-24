# Dirichlet process to draw a random distribution F from a prior F0
#  alpha determines the similarity with the prior guess
#  F0 is the cdf of the prior guess (an object of class "DiscreteDistribution")

##Useful Functions
library(distr) 

cdf_sample <- function(emp_cdf, n=1e3) {
  emp_cdf@r(n) ##random sample
}

dp <- function(alpha, F0, n=1e3) { # n should be large since it's an approx for +oo
  
  s <- cdf_sample(F0,n)            # step 1: draw from F0
  V <- rbeta(n,1,alpha)            # step 2: draw from beta(1,alpha)
  w <- c(1, rep(NA,n-1))           # step 3: compute 'stick breaking process' into w[i]
  w[2:n] <- sapply(2:n, function(i) V[i] * prod(1 - V[1:(i-1)]))
  
  # return the sampled function F which can be itself sampled 
  # this F is a probability mass function where each s[i] has mass w[i]
  function (size=1e4) {
    sample(s, size, prob=w, replace=TRUE)
  }
}

##GSDP Function - MCMC Posterior Simulation Model
gsdp<-function(y,x,
                 spatial=T,varyparam=rep(F,3),nk=5,
                 iters=5000,burn=1000,verbose=10,thin=1,mx.siga,mx.sigb,mx.taua,mx.taub){
  #y:         data (ns x nt)
  #x:         spatial coordinates (ns x 2)
  #spatial:   model the resids as spatially correlated?
  #varyparam: allow the loc, scale and shape to vary by location?
  #nk:        the number of mixture components in the DP (nk=1
  #               is the usual Gaussian copula)z
  #iters:     number of MCMC samples
  #burn:      number of samples to discard
  #verbose:   how often to make updates
  #thin:      Degree of thinning
  ns<- nrow(y) #number of locations
  nt<- ncol(y) #number of years/replications
  xdim<-ncol(x) #covariate dimension
  
  #initial values
  beta<-rep(0,xdim) #initializing beta
  beta<-solve(t(x)%*%x)%*%t(x)%*%rowMeans(y,na.rm=TRUE)
  
  siga<- rep(mx.siga/2, nt)#taking uninformative priors for a's and b's
  sigb<-rep(mx.sigb/2,nt) #TODO: FIGURE THIS OUT
  taua<-rep(mx.taua/2,nt)
  taub<-rep(mx.taub/2,nt)
  bphi<-runif(1,0,.5) #uninformative prior for phi, using .1 based on (Gelfand 2005: 1024)
  phi<-runif(1,0,bphi)
  
  alpha<-rgamma(1,1,1)
  tau<-1/rgamma(1,taua,taub) #initialize tau
  sigma<-1/rgamma(1,siga,sigb) #initialize sigma
  H<-matrix(0,nrow=ns,ncol=ns)#produce the correlation matrix
  for (i in 1:ns){
    for (j in 1:ns){
      H[i,j]<-exp(-phi*(sqrt((x[i,1] - x[j,1])^2 + (x[i,2]- x[j,2])^2)))
    }
  }
  F0<-rmvnorm(1,rep(0,ns),sigma*H)
  theta<-dp(alpha,DiscreteDistribution(F0)) #initialize theta
  
  #MCMC Loop
  
  for (i in 1:iters){
    
    #a)
    #Set Lambda
    Lambda <- solve((tau^(-2) * diag(ns) + sigma^(-2)*solve(H)))
    #Set q0
    q0<-alpha*(det(Lambda))^(1/2)*exp(-0.5*tau^(-2)*(t(y[,t]-x%*%beta)%*%(diag(ns) - tau^(-2)*Lambda)%*%(y[,t]-x%*%beta))) * ((2*pi*tau^2*sigma^2)^(ns/2)*det(H)^(1/2))^(-1)
    theta_t <- rmvnorm(1, tau^(-2)*(Lambda*(rowMeans(y,na.rm=TRUE) - x%*%beta)),Lambda)
  }
}


