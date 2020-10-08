library(MCMCpack)
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
  #s <- sample(1:n.terms,n.terms)  #take draws from F0
  F0.vec <- t(F0)
  theta.unique <- theta.exist
  
  Fnew <- matrix(NA, ncol = ncol(F0.vec), nrow = length(F0.vec[,1]))
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
  
  index <- sample(1:n.terms, size, prob = w, replace=TRUE) #apply stick breaking probabilities to the draws
  Fnew[,index]
}

#compute the stick-breaking weights
makeprobs<-function(v){ 
  N<-length(v)
  probs<-v
  probs[2:N]<-probs[2:N]*cumprod(1-v[2:N-1])
  probs}

#generate categorical variables:
rcat<-function(ndraws,prob){
  prob[is.na(prob)]<-0
  (1:length(prob))%*%rmultinom(ndraws,1,prob)
}


#Function w/ cluster analysis
gsdp.cluster <- function(y,x,
               spatial=T,varyparam=rep(F,3),nk=5,
               iters=1000,burn=400,verbose=10,thin=1,mx.siga,mx.sigb,mx.taua,mx.taub, mx.bphi = 0.1,
               a.alpha = 1, b.alpha = 1,
               display = 5, K = 25) {
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
  #K:         Prior number of spatial clusters
  
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
  
  points.1 <- rep(NA, iters)
  
  ##Priors for spatial clustering (weights and state space model)
  clustervec <- rep(0, ns)
  K <- 5 #prior for finite cluster
  U.cat <- rbeta(K,1,.25) #prior for b.U
  clustervec <- rcat(K, U.cat)
  
  #Prior for number of clusters
  K <- length(unique(clustervec))
  beta.vec <- c(0, rnorm( n = (K - 1), mean = 0, sd = 1)) ##NOT SURE ABOUT THE VARIANCE PRIOR HERE
  G <- matrix(rep(0, K*K), ncol = K)
  g.diag <- runif(K, -1, 1)
  diag(G) <- g.diag
  
  #Priors for phi term
  rho.vec <- c(0,runif(K-1, -1, 1))
  cluster.scale <- 1
  eta.Sigma <- diag(rgamma(K,2,1)^(-1))
  phi.mat <- matrix(rep(0, ns*K), ncol = K)
  lambda.vec <- rgamma(K,2,1)^(-1)
  
  phi.mat[,1] <- rmvnorm(1, rep(0,ns), lambda.vec[k]*H)
  for (k in 2:K) {
    phi.mat[,k] <- rho.vec[k]*phi.mat[,k-1] + rmvnorm(1, rep(0,ns), lambda.vec[k]*H)
  }
  
  #Need prior for the z??
  
  ##MCMC Algorithm
  for (i in 1:iters){
    
    #a)
    
    ##FULL CONDITIONALS
    #Set Lambda
    Lambda <- solve((tau^(-2)*diag(ns) + sigma^(-2)*solve(H)))
    
    #Set q0, h
    q0 <- alpha
    
    htheta<-rmvnorm(10000, tau^(-2)*(Lambda%*%(rowMeans(y) - x%*%t(beta))), Lambda)
    
    #Generate theta vector (DP SAMPLING)
    theta.post <- dp.post(alpha = q0, F0 = htheta, n = sum(T.ast), T.vec = T.ast, theta.exist =theta.exist)[, sample(1:1000,nt)]
    theta.exist <- unique(matrix(c(theta.exist,unique(theta.post,MARGIN = 2)),nrow = ns), MARGIN = 2) #redefine effects
    
    T.length.old <- length(T.ast)
    T.ast<- c(T.ast, rep(0, max(ncol(theta.exist) - T.length.old,0)))
    for (k in 1:ncol(theta.exist)) {
      for (j in 1:ncol(theta.post)) {
        if(mean(theta.exist[,k] == theta.post[,j]) == 1)
          T.ast[k] <- T.ast[k] + 1
      }
    }
    
    #b)

    #c)

    sum.T.ast<-length(T.ast)
    
    Sigma_b <- sd(y)^(2)*solve(t(x)%*%x)
    Sigma_tilde <- solve(solve(Sigma_b)+tau^(-2)*nt*t(x)%*%x)
    sum_d <- 0
    for (j in 1:nt) {
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
    
    ##Spatial Clustering using the State Space Model
    #y_t = H_t z_t + e_t
    #z_t = G z_t-1 + eta_t
    #Posterior Computation
    weights <- matrix(rep(0, ns*K), ncol = K)
    H.mat <- matrix(rep(0,ns*K), ncol = K)
    for (t in 1:nt) {
      weights[t,] <- exp(1*beta.vec + phi.mat[t,])/sum(exp(1*beta.vec + phi.mat[t,]))
      H.mat[t,] <- t(rmultinom(1, 1, prob = weights[t,]))
    }
    #Unsure what G is
    
    
    clustervec.post <- c(1)
    
    #given some cluster count K, we have K different beta's, and cluster specific spatio-temporal random effects (K number of these)
    #beta.vec is some px1 vector at each replicate t
    #x(s) in Paci-Finazzi (2017) is the vector of covariates, here we use y (at each t)
    #phi.vec is the Autoregressive time series, where the random term comes from the independently distributed spatial covariance function
    #K is number of clusters in this iteration
    
    ###MAPPING
    ##quilt.plot()
    ##library(LatticeKrig) # quilt.plot
    ##library(maps) # map()
    
    
    
    
    #rbPal <- colorRampPalette(c('red','blue'))
    theta.pred <- theta.exist[,sample(seq(from = 1, to = length(T.ast), by = 1), size = 1,prob = T.ast/sum(T.ast))]
    #thetaCol <- rbPal(10)[as.numeric(cut(x %*% t(beta) + theta.pred ,breaks = 10))]
    
    print(clustervec.post)
    
 #   count<-count + 1
 #   #hist(y[4,], freq = FALSE, xlim = c(-1,1))
 #   if(count == display){
 #     #par(mfrow = c(2,1))
 #     plot(density(points.1[1:i]))
 #     count <- 0
 #   }
 #   
 # }
 # 
 # plot(density(points.1[burn:iters]))
  
  }
  
}
