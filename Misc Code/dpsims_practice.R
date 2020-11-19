library(mvtnorm)
library(distr)
library(MCMCpack)

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
dp.post <-function(alpha, F0, n, T.vec, theta.exist, n.terms = 1e3,size = 1e4) {
  alpha.post <- alpha + n   #new alpha.post
  s <- sample(1:n.terms,n.terms)  #take draws from F0
  F0.vec <- t(F0[s,])
  theta.unique<- theta.exist

  Fnew <- matrix(NA,ncol = ncol(F0.vec),nrow = length(F0.vec[,1]))
  Fnew[,1] <- F0.vec[,1] #where F0 is the base distribution you are using, probably a normal distribution
  theta.unique <- matrix(c(theta.unique, Fnew[,1]), nrow = length(F0.vec[,1]))
  for(i in 2:ncol(F0.vec)) {
    u <- runif(1)
    if(u < alpha/alpha.post){
      Fnew[,i] <- F0.vec[,sample(1:ncol(F0.vec),size = 1)] #(as before)
      T.vec <- c(T.vec,1)
      theta.unique <- matrix(c(theta.unique, Fnew[,i]),nrow = length(F0.vec[,1]))
    } else {
    # Here, you sample one of the previous locations with probabilities which are proportional to their frequencies; it should not matter if the "probabilities" are larger than 1, because R should renormalises them

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

sigma<-1
alpha<-10

F0<-rmvnorm(10000,rep(0,10),diag(sigma,10))
theta<-dp.mv(alpha,F0)

theta.sample<-theta[,sample(1:1000,100)]


theta.exist <- unique(theta.sample,MARGIN = 2)
T.ast<-rep(0, ncol(theta.exist))
for (i in 1:ncol(theta.exist)){
  for (j in 1:ncol(theta.sample)) {
    if(theta.exist[,i] == theta.sample[,j])
       T.ast[i] <- T.ast[i] + 1
  }
}

theta.post <- dp.post(alpha, F0, n = sum(T.ast), T.vec = T.ast, theta.exist = theta.exist)

par(mfrow = c(3,1))

hist(colMeans(theta))
hist(colMeans(theta.sample))
hist(colMeans(theta.post))

###ORIGINAL DP FUNCTION

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
