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

##function without cluster implementation
gsdp_no_cluster<-function(y,x,
                 spatial=T,varyparam=rep(F,3),nk=5,
                 iters=1000,burn=100,verbose=10,thin=1,mx.siga,mx.sigb,mx.taua,mx.taub, noClusters = NA){
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
  #noClusters: number of clusters
  
  ns<- nrow(y) #number of locations
  nt<- ncol(y) #number of years/replications
  xdim<-ncol(x) #covariate dimension
  noClusters<-ns
  
  #initial values
  beta<-rep(0,xdim) #initializing beta
  beta<-solve(t(x)%*%x)%*%t(x)%*%rowMeans(y,na.rm=TRUE)
  
  siga<- rep(mx.siga/2, nt)#taking uninformative priors for a's and b's
  sigb<-rep(mx.sigb/2,nt) #TODO: FIGURE THIS OUT
  taua<-rep(mx.taua/2,nt)
  taub<-rep(mx.taub/2,nt)
  bphi<-runif(1,0,.5) #uninformative prior for phi, using .1 based on (Gelfand 2005: 1024)
  phi<-runif(1,0,bphi)
  
  alpha<-rgamma(1,1,1) #a_alpha and b_alpha set as 1's
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
  
  htheta<-matrix(0,nrow=ns,ncol=nt)
  wvec<-rep(0,nt)
  q0<-rep(0,nt)
  qj<-matrix(0,nrow=ns,ncol=nt)
  denomvec<-rep(0,nt)
  thetavec<-matrix(0,nrow=ns,ncol=nt)
  thetaspec<-matrix(0,nrow=ns,ncol=nt)
  
  for (i in 1:iters){
    
    #a)
    #Set Lambda
    Lambda <- solve((tau^(-2) * diag(ns) + sigma^(-2)*solve(H)))
    
    #Set q0, h
    for (j in 1:nt){
      q0[j]<-alpha*(det(Lambda))^(1/2)*exp(-0.5*tau^(-2)*(t(y[,j]-x%*%t(beta))%*%(diag(ns) - tau^(-2)*Lambda)%*%(y[,j]-x%*%t(beta)))) * ((2*pi*tau^2*sigma^2)^(ns/2)*det(H)^(1/2))^(-1)
      htheta[,j] <- q0[j] * rmvnorm(1, tau^(-2)*(Lambda%*%(y[,j] - x%*%t(beta))),Lambda)
    }
    for (j in 1:nt){
      qj[,j]<-rmvnorm(1, x%*%t(beta) + htheta[,j],tau^(2)*diag(ns))
      #NOTE: we use htheta here we we do not introduce clustering in this code, s.t. each location is 'distinct'
    }
    for(j in 1:nt){
      denomvec[j]<-q0[j]+sum(rowSums(qj))
    }
    
    #posterior theta
    for (j in 1:nt){
      thetavec[,j]<-  (q0[j]*htheta[,j] + qj[,j])/(denomvec[j])
    }
    
    #b)
    err<-matrix(0,nrow=ns,ncol=nt)
    for (j in 1:nt){
      err[,j]<-y[,j] - x%*%t(beta)
    }
    for (j in 1:noClusters){
      thetaspec[,j]<-rmvnorm(1,tau^(-2)*solve(1*tau^(-2)*diag(ns)+sigma^(-2)*solve(H))%*% err[,j],solve(1*tau^(-2)*diag(ns)+sigma^(-2)*solve(H)))
      #NOTE: 1* is normally T_j, not implementing clustering
    }
    
    #c)
    Sigma_b <-sd(y)^(2)*solve(t(x)%*%x)
    Sigma_tilde<-solve(solve(Sigma_b)+tau^(-2)*nt*t(x)%*%x)
    #NOTE *nt* is here instead of sum of T, using stationary covariates
    sum_d <-0
    for (j in 1:nt){
      sum_d<-sum_d + t(x)%*%(y[,j] - thetavec[,j])
    }
    beta_tilde<-Sigma_tilde%*%(solve(Sigma_b)%*%t(beta) + tau^(-2)*sum_d)
    beta<-rmvnorm(1,beta_tilde,Sigma_tilde)
    
    sum_d2<-0
    for (j in 1:nt){
      sum_d2 <- sum_d2 + t(y[,j] - x%*%t(beta) - thetavec[,j])%*%(y[,j] -x%*%t(beta) - thetavec[,j])
    }
    tau<-rgamma(1,taua+0.5*nt*noClusters,taub+0.5*sum_d2[1])
    
    #d)
    #update alpha
    eta_aux<-rbeta(1,alpha+1,nt)
    p_d<-(1+noClusters-1)/(nt*(1-log(eta_aux))+1+noClusters -1) #a_alpha and b_alpha set as 1
    alpha<-p_d*rgamma(1,1+noClusters,1-log(eta_aux)) + (1-p_d)*rgamma(1,1+noClusters-1,1-log(eta_aux))
    #update sigma
    sum_d3<-0
    for (j in 1:nt){
      sum_d3 <- sum_d3 + t(thetaspec[,j])%*%solve(H)%*%thetaspec[,j]
    }
    sigma<-rgamma(1,siga+.5*nt*noClusters,sigb+.5*sum_d3[1])
    #update phi
    phi<-(det(H))^(-noClusters/2)*exp(-sum_d3/(2*sigma^(2)))
  }
}


