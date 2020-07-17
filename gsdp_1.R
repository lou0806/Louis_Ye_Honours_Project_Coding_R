# Dirichlet process to draw a random distribution F from a prior F0
#  alpha determines the similarity with the prior guess
#  F0 is the cdf of the prior guess (an object of class "DiscreteDistribution")

##Useful Functions
library(distr)
library(mvtnorm)

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
##This function uses in-model Kriging (discussed in end of Section 3 in Gelfand (2005)) to 'fill in' missing data

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

  #initial values
  beta<-rep(0,xdim) #initializing beta
  beta<-solve(t(x)%*%x)%*%t(x)%*%rowMeans(y,na.rm=TRUE)
  beta<-t(beta)
  
  siga<- mx.siga/2#rep(mx.siga/2, nt)#taking uninformative priors for a's and b's
  sigb<-mx.sigb/2#rep(mx.sigb/2,nt) #TODO: FIGURE THIS OUT
  taua<-mx.taua/2#rep(mx.taua/2,nt)
  taub<-mx.taub/2#rep(mx.taub/2,nt)
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
  
  #MCMC Loop initialize
  htheta<-matrix(0,nrow=ns,ncol=nt)
  q0<-rep(0,nt)
  qj<-matrix(0,nrow=ns,ncol=nt)
  denomvec<-rep(0,nt)
  
  wvec<-matrix(NA,nrow=ns,ncol=nt)
  thetavec<-matrix(0,nrow=ns,ncol=nt)
  clustercount <- rep(NA, nt)
  thetaspec <- matrix(NA,nrow=ns,ncol=nt)
  groupvec <- matrix(NA,nrow=ns,ncol=nt)
  
  #Within-replicate Kriging using the initialized theta, tau and beta
  if(anyNA(y) == TRUE) {
    for (i in 1:nt) {
      x.u<-matrix(c(x[,1][is.na(y[,i])],x[,2][is.na(y[,i])]),ncol=2)
      y.obs<-y[,i][!is.na(y[,i])]
      x.obs<-matrix(c(x[,1][!is.na(y[,i])],x[,2][!is.na(y[,i])]),ncol=2)
      n.na<-length(y[,i]) - length(y.obs)
      beta.temp<-rep(0,length(y.obs)) #initializing beta
      beta.temp<-solve(t(x.obs)%*%x.obs)%*%t(x.obs)%*%y.obs
      #beta.temp<-t(beta.temp)
      thetavec.temp <- sample(theta(),n.na)
      y.u<-rmvnorm(1,x.u%*%beta.temp + thetavec.temp,tau^2*diag(n.na))
      
      y[,i][is.na(y[,i])]<-y.u
    }
  }
  
  
  ##MCMC Algorithm
  for (i in 1:iters){
    
    #a)
    #Generate theta vector
    cluster_count<-0
    for (i in 1:nt) {
      thetavec[,i] <- sample(theta(),ns)
      groupvec[,i] <- c(unique(thetavec[,i]),rep(NA, ns - length(unique(thetavec[,i]))))
      cluster_count[i] <- sum(groupvec[,i] > 0,na.rm = TRUE)
      for (j in 1:cluster_count[i]) {
        wvec[thetavec == groupvec[j]] <- j
      }
    }
    theta.ast <- unique.array(thetavec, MARGIN = 2)
    T.ast <- length(unique.array(thetavec, MARGIN = 2))/ns
    
    ##FULL CONDITIONALS - below code can be useful later
    #Set Lambda
#    Lambda <- solve((tau^(-2)*diag(ns) + sigma^(-2)*solve(H)))
    
    #Set q0, h
#    for (j in 1:nt){
#      #q0[j]<-alpha*(det(Lambda))^(1/2) * exp(-0.5*tau^(-2)*(t(y[,j]-x%*%t(beta))%*%(diag(ns) - tau^(-2)*Lambda)%*%(y[,j]#-x%*%t(beta)))) * ((2*pi*tau^2*sigma^2)^(ns/2)*det(chol(H)))^(-1)
#      #Using log
#      q0.temp <- log(alpha*det(chol(Lambda))) + -0.5*tau^(-2)*(t(y[,j]-x%*%t(beta))%*%(diag(ns) - tau^(-2)*Lambda)%*%(y[#,j]-x%*%t(beta))) + log(((2*pi*tau^2*sigma^2)^(ns/2)*det(chol(H)))^(-1))
#      q0[j] <- exp(q0.temp)
#      htheta[,j] <- q0[j] * rmvnorm(1, tau^(-2)*(Lambda%*%(y[,j] - x%*%t(beta))),Lambda)
#      
#      qj[,j]<-rmvnorm(1, x%*%t(beta) + htheta[,j],tau^(2)*diag(ns))
#      denomvec[j]<-q0[j]+sum(rowSums(qj))
#      thetavec[,j]<-  (1/(denomvec[j]))*(q0[j]*htheta[,j] + qj[,j])
#    }
      
    ##FULL CONDITIONALS
#    err<-matrix(0,nrow=ns,ncol=nt)
#    for (j in 1:nt){
#      err[,j]<-y[,j] - x%*%t(beta)
#    }
#    for (j in 1:noClusters){
#      thetaspec[,j]<-rmvnorm(1,tau^(-2)*solve(1*tau^(-2)*diag(ns)+sigma^(-2)*solve(H))%*% err[,j],solve(1*tau^(-2)*diag#(ns)+sigma^(-2)*solve(H)))
#    }
    
    #c)
    ##TODO: IMPLEMENT THE CLUSTERING and GENERATE NEW SURFACE (Gelfand (2005) pg 1025)
    newTheta <- thetavec[,ceiling(runif(1,0,nt))] ##UNSURE WHERE TO SAMPLE FROM
    sum.T.ast <- 0
    for (i in 1:nt) {
      if (mean(thetavec[,i] == newTheta) == 1) {
        sum.T.ast <- sum.T.ast + 1
      }
    }
    theta_new <- alpha/(alpha + nt) * F0 + 1/(alpha + nt) * sum.T.ast
    #clust_new <- unique(theta_new)
    
    Sigma_b <-sd(y)^(2)*solve(t(x)%*%x)
    Sigma_tilde<-solve(solve(Sigma_b)+tau^(-2)*nt*t(x)%*%x)
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
    tau<-1/rgamma(1,taua+0.5*nt*ns,taub+0.5*sum_d2[1])
    
    #d)
    #update alpha
    eta_aux<-rbeta(1,alpha+1,nt)
    p_d<-(1+sum.T.ast-1)/(nt*(1-log(eta_aux))+1+sum.T.ast -1) #a_alpha and b_alpha set as 1
    alpha<-p_d*rgamma(1,1+sum.T.ast,1-log(eta_aux)) + (1-p_d)*rgamma(1,1+sum.T.ast-1,1-log(eta_aux))
    #update sigma
    sum_d3<-0
    for (j in 1:T.ast){
      sum_d3 <- sum_d3 + t(theta.ast[,j])%*%solve(H)%*%theta.ast[,j]
    }
    sigma<-rgamma(1,siga+.5*nt*sum.T.ast,sigb+.5*sum_d3)
    #update phi
    phi <- runif(1,0,bphi)
    #phi<-(det(chol(H)))^(-sum.T.ast)*exp(-sum_d3/2*sigma^2)
  
    #for (i in 1:ns){
    #  for (j in 1:ns){
    #    H[i,j]<-exp(-phi*(sqrt((x[i,1] - x[j,1])^2 + (x[i,2]- x[j,2])^2)))
    #  }
    #}
    F0<-rmvnorm(1,rep(0,ns),sigma*H)
    theta<-dp(alpha,DiscreteDistribution(F0)) #redo theta
  }
}


