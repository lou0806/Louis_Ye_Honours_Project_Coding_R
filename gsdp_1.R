# Dirichlet process to draw a random distribution F from a prior F0
#  alpha determines the similarity with the prior guess
#  F0 is the cdf of the prior guess (an object of class "DiscreteDistribution")

##Useful Functions
library(distr)
library(mvtnorm)


#Custom DP
dp.mv <- function(alpha, F0, n=1e3, size = 1e4) { # n should be large since it's an approx for +oo
  
  s <- sample(1:n,n)            # step 1: take draws from F0
  F0.vec <- t(F0[s,])
  
  V <- rbeta(n,1,alpha)            # step 2: draw from beta(1,alpha)
  
  w <- c(1, rep(NA,n-1))           # step 3: compute 'stick breaking process' into w[i]
  w[2:n] <- sapply(2:n, function(i) V[i] * prod(1 - V[1:(i-1)]))
  
  index <- sample(s, size, prob= w, replace=TRUE) #apply stick breaking probabilities to the draws
  F0.vec[,index]
}

#Posterior sampling
#T.vec represents the vector of cluster sizes
#Weighted average interpretation
dp.post <-function(alpha, F0, n, T.vec, theta.exist, n.terms = 1e3, size = 1e4) {
  alpha.post <- alpha + n   #new alpha.post
  s <- sample(1:n.terms,n.terms)  #take draws from F0
  F0.vec <- t(F0[s,])
  theta.unique<- theta.exist
  
  Fnew <- matrix(NA,ncol = ncol(F0.vec),nrow = length(F0.vec[,1]))
  Fnew[,1] <- F0.vec[,1]
  theta.unique <- matrix(c(theta.unique, Fnew[,1]), nrow = length(F0.vec[,1]))
  for(i in 2:ncol(F0.vec)) {
    u <- runif(1)
    if(u < alpha/alpha.post){
      Fnew[,i] <- F0.vec[,sample(1:ncol(F0.vec),size = 1)] #(as before)
      T.vec <- c(T.vec,1)
      theta.unique <- matrix(c(theta.unique, Fnew[,i]),nrow = length(F0.vec[,1]))
    } else {
      
      index.temp <- sample(1:length(T.vec),size = 1, prob = c(T.vec))
      T.vec[index.temp] <-T.vec[index.temp] + 1
      Fnew[,i] <- theta.unique[,index.temp]
    }
  }
  
  V <- rbeta(n.terms, 1,alpha.post) #draw from beta(1,alpha.post)
  w <- c(1, rep(NA,n.terms-1))
  w[2:n.terms] <- sapply(2:n.terms, function(i) V[i] * prod(1 - V[1:(i-1)]))
  
  index <- sample(s, size, prob= w, replace=TRUE) #apply stick breaking probabilities to the draws
  Fnew[,index]
}

##GSDP Function - MCMC Posterior Simulation Model
##This function uses in-model Kriging (discussed in end of Section 3 in Gelfand (2005)) to 'fill in' missing data

##function without cluster implementation
gsdp<-function(y,x,
                 spatial=T,varyparam=rep(F,3),nk=5,
                 iters=1000,burn=50,verbose=10,thin=1,mx.siga,mx.sigb,mx.taua,mx.taub, mx.bphi = 0.1,
                 a.alpha = 1, b.alpha = 1,
                 noClusters = NA,display = 5){
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
  
  count<-0
  clustercount<-rep(0,iters)

  #initial values
  beta<-rep(0,xdim) #initializing beta
  beta<-solve(t(x)%*%x)%*%t(x)%*%rowMeans(y,na.rm=TRUE)
  beta<-t(beta)
  
  siga<- mx.siga/2#rep(mx.siga/2, nt)#taking uninformative priors for a's and b's
  sigb<-mx.sigb/2#rep(mx.sigb/2,nt) #TODO: FIGURE THIS OUT
  taua<-mx.taua/2#rep(mx.taua/2,nt)
  taub<-mx.taub/2#rep(mx.taub/2,nt)
  bphi<-runif(1,0,mx.bphi) #uninformative prior for phi, using .1 based on (Gelfand 2005: 1024)
  phi<-runif(1,0,bphi)
  
  alpha<-rgamma(1,a.alpha,b.alpha) #a_alpha and b_alpha set as 1's
  ##Create vector to save alpha across iterations
  
  tau<-1/rgamma(1,taua,taub) #initialize tau
  sigma<-1/rgamma(1,siga,sigb) #initialize sigma
  H<-matrix(0,nrow=ns,ncol=ns)#produce the correlation matrix
  for (i in 1:ns){
    for (j in 1:ns){
      H[i,j]<-exp(-phi*(sqrt((x[i,1] - x[j,1])^2 + (x[i,2]- x[j,2])^2)))
    }
  }
  F0<-rmvnorm(10000,rep(0,ns),sigma*H)
  theta<-dp.mv(alpha = alpha, F0 = F0)[,sample(1:1000,nt)] #initialize theta, ~ iid G^(n)
  theta.exist <- unique(theta,MARGIN = 2) #initialize clusters
  T.ast<-rep(0, ncol(theta.exist))
  for (i in 1:ncol(theta.exist)){
    for (j in 1:ncol(theta)) {
      if(mean(theta.exist[,i] == theta[,j]) == 1)
        T.ast[i] <- T.ast[i] + 1
    }
  }
  theta.post<-theta.exist[,1]
  
  #MCMC Loop initialize
  htheta<-matrix(0,nrow=ns,ncol=nt)
  q0<-rep(0,nt)
  qj<-matrix(0,nrow=ns,ncol=nt)
  denomvec<-rep(0,nt)
  
  wvec<-matrix(NA,nrow=ns,ncol=nt)
  thetavec<-matrix(0,nrow=ns,ncol=nt)
  thetaspec <- matrix(NA,nrow=ns,ncol=nt)
  groupvec <- matrix(NA,nrow=ns,ncol=nt)
  
  #Within-replicate Kriging using the initialized theta, tau and beta
  #if(anyNA(y) == TRUE) {
  #  for (i in 1:nt) {
  #    x.u<-matrix(c(x[,1][is.na(y[,i])],x[,2][is.na(y[,i])]),ncol=2)
  #    y.obs<-y[,i][!is.na(y[,i])]
  #    x.obs<-matrix(c(x[,1][!is.na(y[,i])],x[,2][!is.na(y[,i])]),ncol=2)
  #    n.na<-length(y[,i]) - length(y.obs)
  #    beta.temp<-rep(0,length(y.obs)) #initializing beta
  #    beta.temp<-solve(t(x.obs)%*%x.obs)%*%t(x.obs)%*%y.obs
  #    #beta.temp<-t(beta.temp)
  #    thetavec.temp <- sample(theta(),n.na)
  #    y.u<-rmvnorm(1,x.u%*%beta.temp + thetavec.temp,tau^2*diag(n.na))
  #    
  #    y[,i][is.na(y[,i])]<-y.u
  #  }
  #}
  
  
  ##MCMC Algorithm
  for (i in 1:iters){
    
    #a)
    
    ##FULL CONDITIONALS
    #Set Lambda
    Lambda <- solve((tau^(-2)*diag(ns) + sigma^(-2)*solve(H)))
    
    #Set q0, h
#    q0<-alpha*(det(Lambda))^(1/2) * exp(-0.5*tau^(-2)*(t(y[,j]-x%*%t(beta))%*%(diag(ns) - tau^(-2)*Lambda)%*%(y[,j]-x%*%t(beta)))) * ((2*pi*tau^2*sigma^2)^(ns/2)*det(chol(H)))^(-1)
#    q0.temp <- log(alpha*det(chol(Lambda))) + -0.5*tau^(-2)*(t(y[,j]-x%*%t(beta))%*%(diag(ns) - tau^(-2)*Lambda)%*%(y[,j]-x%*%t(beta))) + log(((2*pi*tau^2*sigma^2)^(ns/2)*det(chol(H)))^(-1))
    q0 <- alpha
    
    htheta<-rmvnorm(10000, tau^(-2)*(Lambda%*%(rowMeans(y) - x%*%t(beta))), Lambda)
    
    #Generate theta vector (DP SAMPLING)
    ## TODO: USE DPpackage, likely LDDPdensity()
   #cluster_count<-0
   #for (k in 1:nt) {
   #  thetavec[,k] <- sample(theta(),ns)
   #  groupvec[,k] <- c(unique(thetavec[,k]),rep(NA, ns - length(unique(thetavec[,k]))))
   #  cluster_count[k] <- sum(groupvec[,k] > 0,na.rm = TRUE)
   #  for (j in 1:cluster_count[k]) {
   #    wvec[thetavec == groupvec[j]] <- j
   #  }
   #}
    
    ##TODO: Use weigths h(theta) and q0 in the posterior DP sampling
    theta.post <- dp.post(q0,htheta,sum(T.ast),n = sum(T.ast),T.vec = T.ast, theta.exist =theta.exist)[,sample(1:1000,nt)]
    theta.exist <- unique(matrix(c(theta.exist,unique(theta.post,MARGIN = 2)),nrow = ns), MARGIN = 2) #redefine effects
    T.length.old <- length(T.ast)
    T.ast<- c(T.ast, rep(0, max(ncol(theta.exist) - T.length.old,0)))
    for (i in 1:ncol(theta.exist)){
      for (j in 1:ncol(theta.post)) {
        if(mean(theta.exist[,i] == theta.post[,j]) == 1)
          T.ast[i] <- T.ast[i] + 1
      }
    }
    
    
    
    #b)
    #theta.ast <- unique.array(thetavec, MARGIN = 2)
    #T.ast <- length(unique.array(thetavec, MARGIN = 2))/ns
  
    #c)
    ##TODO: IMPLEMENT THE CLUSTERING and GENERATE NEW SURFACE (Gelfand (2005) pg 1025)
    #newTheta <- thetavec[,ceiling(runif(1,0,nt))] ##UNSURE WHERE TO SAMPLE FROM
    #sum.T.ast <- 0
    #for (j in 1:nt) {
    #  if (mean(thetavec[,j] == newTheta) == 1) {
    #    sum.T.ast <- sum.T.ast + 1
    #  }
    #}
    #theta_new <- alpha/(alpha + nt) * F0 + 1/(alpha + nt) * sum.T.ast
    #clust_new <- unique(theta_new)
    
    sum.T.ast<-length(T.ast)
    
    Sigma_b <-sd(y)^(2)*solve(t(x)%*%x)
    Sigma_tilde<-solve(solve(Sigma_b)+tau^(-2)*nt*t(x)%*%x)
    sum_d <-0
    for (j in 1:nt){
      sum_d<-sum_d + t(x)%*%(y[,j] - theta.post[,j])
    }
    
    beta_tilde<-Sigma_tilde%*%(solve(Sigma_b)%*%t(beta) + tau^(-2)*sum_d)
    beta<-rmvnorm(1,beta_tilde,Sigma_tilde)
    
    sum_d2<-0
    for (j in 1:nt){
      sum_d2 <- sum_d2 + t(y[,j] - x%*%t(beta) - theta.post[,j])%*%(y[,j] -x%*%t(beta) - thetavec[,j])
    }
    tau<-1/rgamma(1,taua+0.5*nt*ns,taub+0.5*sum_d2[1])
    
    #d)
    #update alpha
    eta_aux<-rbeta(1,alpha+1,nt)
    p_d<-(1+sum.T.ast-1)/(nt*(1-log(eta_aux))+1+sum.T.ast -1) #a_alpha and b_alpha set as 1
    alpha<-p_d*rgamma(1,1+sum.T.ast,1-log(eta_aux)) + (1-p_d)*rgamma(1,1+sum.T.ast-1,1-log(eta_aux))
    
    #update sigma
    sum_d3<-0
    for (j in 1:length(T.ast)){
      sum_d3 <- sum_d3 + t(theta.exist[,j])%*%solve(H)%*%theta.exist[,j]
    }
    sigma<-rgamma(1,siga+.5*nt*sum.T.ast,sigb+.5*sum_d3)
    #update phi and H
    phi <- runif(1,0,bphi)
    for (k in 1:ns){
      for (j in 1:ns){
        H[k,j]<-exp(-phi*(sqrt((x[k,1] - x[j,1])^2 + (x[k,2]- x[j,2])^2)))
      }
    }
    theta<-theta.post #redo theta
    
    #Iteration plots
    #gvec<-rep(NA,ns)
    #for (j in 1:ns) {
    #  gvec[j] = which(newTheta[j] == unique(newTheta))
    #}
    
    
    ##
    #clustercount[i] <- length(T.ast)
    
    #DEBUG
    #if(clustercount[i] == 0) {
    #  print(newTheta)
    #}
    count<-count + 1
    
    if(count==display){
      par(mfrow = c(2,1))
      #plot(clustercount)
      plot(T.ast > 1)
      plot(T.ast)
      count <- 0
    }
    
  }
  print(theta.unique)
}


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


#old phi update

#phi<-(det(chol(H)))^(-sum.T.ast)*exp(-sum_d3/2*sigma^2)

#for (i in 1:ns){
#  for (j in 1:ns){
#    H[i,j]<-exp(-phi*(sqrt((x[i,1] - x[j,1])^2 + (x[i,2]- x[j,2])^2)))
#  }
#}
