#############################################################################
# Do MCMC sampling for the spatial stick-breaking model
#    Here's the model:
#    
#    y[s] ~ N(x*beta+mu[g[s]],sige[g[s])   
#    mu[k]~dnorm(sigs), k=1,...n.terms
#    sige[1],...,sige[n.terms],sigs~dunif(0,mx.sig)
#
#    g[s] ~ Categorical(p_1[s],...,p_m[s])
#    p_j[s] = v_j[s] * prod_{k<j} (1-v_k[s])
#    v_k[s] = w_j(knots[j,],x[s,],rho) * v[k]
#            (w is a Gaussian kernel)
#
#    knot[j1,j2] ~ unif(0,1)
#    rho ~ Unif(0,1)
#    v[j] ~ beta(1,DOF)
#
#############################################################################

library(MCMCpack)
library(distr)
library(mvtnorm)

#Compute v_k(s)
onevs<-function(z,rho,knot,v,kernel="gaussian"){
  kkk<-v
  for(j in 1:length(v)){
    if(kernel=="gaussian"){
      kkk<-kkk*rho^((z[,j]-knot[j])^2)
    }     
    if(kernel=="uniform"){
      kkk<-kkk*ifelse(abs(z[,j]-knot[j])<rho,1,0.01)
    }     
    
  }
  kkk}

#take all the vs and compute probabilities
makeprobs<-function(vs){
  m<-ncol(vs)
  failures<-matrix(0,m,m)
  for(j in 2:m){
    failures[1:(j-1),j]<-1
  }
  probs<-exp(log(vs)+log(1-vs)%*%failures)
  probs[,m]<-1-apply(probs[,1:m],1,sum)
  probs
}

#generate categorical variables:
rcat<-function(prob){
  (1:length(prob))%*%rmultinom(1,1,ifelse(prob>0,prob,0))
}

exp.H <- function (decay,z) {
  ns <- nrow(z)
  Exp.Mat<-matrix(0,nrow=ns,ncol=ns)#produce the spatial correlation matrix
  
  for (i in 1:ns){
    for (j in 1:ns){
      Exp.Mat[i,j]<-exp(-decay*(sqrt((z[i,1] - z[j,1])^2 + (z[i,2]- z[j,2])^2)))
    }
  }
  return(Exp.Mat)
}

#The whole Shebang!
SSB.dynamic<-function(y,x=NA,z,DOF=1,mx.sige=1,mx.sigs=1,n.terms=100,
                      runs=1000,burn=300,display=10) {
  #y:       data
  #x:       covariates
  #z:       n x 2 matrix of coordinates, scaled to [0,1]
  #n.terms: number of terms in the mixture dist.
  #runs:    number of MCMC samples
  #burn:    number of samples to discard
  #display: how often to display results
  #DOF:     v~beta(1,DOF)
  #mx.sig:  the sds are U(0,mx.sd)
  
  ns <- nrow(y) #number of locations
  nt <- ncol(y) #number of years/replications
  #n.terms <- ncol(y)
  
  n<-ns
  if(is.na(max(x))){x<-matrix(1,n,1)}
  p<-ncol(x)
  
  #Standardize the outcomes and predictors:
  if(min(z)<0 | max(z)>1){print("THE SPATIAL COORDINATES ARE NOT STANDARDIZED!")}   
  
  #initial values
  beta.mat<-matrix(rep(0,nt*p), ncol = nt)
  cluster.count <- rep(0, runs)
  v <- rep(.9,n.terms) #originally .9, can change to influence initial clustering
  sige <- matrix(rep(mx.sige/2,n.terms*nt), nrow = nt)
  taue <- 1/sige^2
  mu.mat <- matrix(rep(0,n.terms*nt), nrow = nt)
  sigs <- rep(mx.sigs/2, nt)
  taumu <- 1/sigs^2
  knot <- matrix(runif(2*n.terms,0,1),2,n.terms)
  rho <- .5 #influences clustering
  
  g.mat <- matrix(rep(1,n*nt), ncol = nt) #membership
  y.mat <- y
  
  vs<-matrix(0,n,n.terms)
  for(k in 1:n.terms){
    vs[,k]<-onevs(z,rho,knot[,k],v[k])
  }
  probs <- makeprobs(vs)
  
  sumtp <- summu <- summu2 <- rep(0,n)
  count <- afterburn <- 0
  keeprho <- rep(0,runs)
  keepbeta <- matrix(0,runs,p)
  knot.mat <- array(rep(knot, nt), dim = c(2, ns, nt))
  
  ##Paci Finazzi terms
  rho.ts <- runif(n = n.terms, min = -1, max = 1) #rho is length of # of clusters
  phi.ts <- matrix(rep(0, nt*n.terms), ncol = nt)
  tau.ts <- rep(0, n.terms)
  
  H<-matrix(0,nrow=ns,ncol=ns)#produce the spatial correlation matrix
  lambda <- rgamma(nt, shape = 1, rate = 2)^{-1}
  decay <- 0.4 #MCMC doesn't update this. Is theta(decay) in Paci-Finazzi
  for (i in 1:ns){
    for (j in 1:ns){
      H[i,j]<-exp(-decay*(sqrt((z[i,1] - z[j,1])^2 + (z[i,2]- z[j,2])^2)))
    }
  }
  
  phi.ts[,1] <- rmvnorm(1,rep(0,ns), rho.ts[1]^2 * H)
  
  #Result terms
  probs.cluster1 <- g.mat#probability of being in cluster 1
  
  for(i in 1:runs) {
    
    #Update beta
    COV <- rep(NA, nt)
    r <- matrix(rep(0, ns*nt), ncol = nt)
    
    for (t in 1:nt) {
      COV[t]<-solve(t(x)%*%diag(taue[t,][g.mat[,t]])%*%x)
      mn<-COV[t]%*%t(x)%*%diag(taue[g.mat[,t]])%*%(y[,t]-mu.mat[g.mat[,t]])
      beta.mat[,t] <- mean(y[,t])
        #mn + t(chol(COV[t]))
        #Altered for no covariate data
      r[,t]<-y[,t]-x%*%beta.mat[,t]
      knot <- knot.mat[,,t]
    
      for(j in 1:n.terms) {
        #update mu
        nobs<-sum(g.mat[,t]==j)
        mu.mat[t,j]<-rnorm(1,0,sigs[t])
        if(nobs>0){
          mu.mat[t,j]<-rnorm(1,mean(r[,t][g.mat[,t]==j]),1/sqrt(taumu[t]+nobs*taue[j]))
        }
        
        #update sige
        cansige<-sige[t,]
        cansige[j] <- rnorm(1,sige[t,j],.1)
        if(sum(g.mat[,t]==j) > 0 & cansige[j]>0 & cansige[j] < mx.sige) {
          MHrate<-sum(dnorm(r[g.mat[,t]==j],mu.mat[t,j],cansige[j],log=T)-
                        dnorm(r[g.mat[,t]==j],mu.mat[t,j],sige[t,j],log=T))
          if(runif(1,0,1) < exp(MHrate)) {
            sige[t,]<-cansige
          } 
        }
        if(sum(g.mat[,t]==j)==0 & cansige[j]>0 & cansige[j]<mx.sige){
          sige[t,]<-cansige 
        }
        taue[t,]<-1/sige[t,]^2
      }

      #update sigs:
      cansigs<-rnorm(1,sigs[t],.1)
      if(cansigs>0 & cansigs<mx.sigs){
        MHrate<-sum(dnorm(mu.mat[t,],0,cansigs,log=T)-
                      dnorm(mu.mat[t,],0,sigs[t],log=T))
        if(runif(1,0,1)<exp(MHrate)){sigs[t]<-cansigs} 
      }
      taumu[t]<-1/sigs[t]^2

      #update g
      for(s in 1:n){
        cang<-g.mat[,t]
        cang[s]<-rcat(probs[s,])
        MHrate<-dnorm(r[s],mu.mat[t,][cang[s]],sige[cang[s]],log=T)-
          dnorm(r[s],mu.mat[t,][g.mat[s,t]],sige[g.mat[s,t]],log=T)
        if(runif(1,0,1)<exp(MHrate)){
          g.mat[,t]<-cang
        }
      }
    
      #update rho:
      canrho<-runif(1,0,1)
      canvs<-vs
      for(k in 1:n.terms){
        canvs[,k]<-onevs(z,canrho,knot[,k],v[k])
      }
      canprobs<-makeprobs(canvs)
      MHrate<-sum(log(canprobs[cbind(1:n, g.mat[,t])])-
                    log(probs[cbind(1:n, g.mat[,t])]))
      if(runif(1,0,1) < exp(MHrate)) {
        rho<-canrho
        vs<-canvs
        probs<-canprobs
      }
    
      #update v:
      for(k in 1:(n.terms-1)) {
        if(max(g.mat[,t])<k) {
          v[k]~rbeta(1,1,DOF)
        }
        if(max(g.mat[,t])>=k){
          canv<-v;canv[k]<-rnorm(1,v[k],.05)
          if(canv[k]>0 & canv[k]<1){
            canvs[,k]<-onevs(z,rho,knot[,k],canv[k])
            canprobs<-makeprobs(canvs)
            MHrate<-sum(log(canprobs[cbind(1:n,g.mat[,t])])-
                          log(probs[cbind(1:n,g.mat[,t])]))+
              log(dbeta(canv[k],1,DOF)/dbeta(v[k],1,DOF))
            if(runif(1,0,1)<exp(MHrate)){
              v<-canv
              vs<-canvs
              probs<-canprobs
            }
          }
        }
      }
      
      #update knots:
      for(j in 1:2) {
        for(k in 1:n.terms) {
          canvs <- vs
          canknot <- knot
          canknot[j,k] <- rnorm(1,knot[j,k],.1)
          if(canknot[j,k] > 0 & canknot[j,k]<1) {
            canvs[,k] <- onevs(z,rho,canknot[,k],v[k])
            canprobs <- makeprobs(canvs)
            MHrate<-sum(log(canprobs[cbind(1:n,g.mat[,t])]) -
                          log(probs[cbind(1:n,g.mat[,t])]))
            if(runif(1,0,1)<exp(MHrate)) {
              knot<-canknot
              vs<-canvs
              probs<-canprobs
            }
          }
        }
      }
      knot.mat[,,t] <- knot
    
      probs.cluster1[,t] <- probs[,1]
    }
    
    #update Paci-Finazzi:
    #separate the components
    #Keeping a track of knots (dynamic)
    n.clust <- rep(0, nt)
    knot.clust <- array(rep(knot, nt), dim = c(2, max(g.mat), nt))
    max.clust <- max(max(g.mat),2)
    phi.ts <- matrix(rep(0, max.clust*nt), ncol = nt)
    phi.effect <- matrix(rep(0, nt*ns), ncol = nt)
    if(length(rho.ts) < max.clust) rho.ts <- c(rho.ts, rep(0, (max.clust - length(rho.ts))))
    
    for (t in 1:nt) {
      n.clust[t] <- length(unique(g.mat[,t]))
      temp.length <- max(g.mat) - ncol(knot.mat[,,t][,unique(g.mat[,t])])
      temp.knot <- cbind(knot.mat[,,t][,unique(g.mat[,t])] , matrix(rep(0.5, (temp.length*2)), nrow = 2))
      knot.clust[,,t] <- temp.knot
    }

    phi.ts[,1] <- rmvnorm(1, rep(0, max.clust), lambda[1]*exp.H(decay, t(knot.clust[,,1] )))
    
    #Unique random effects
    temp.H <- exp.H(decay, t(knot.clust[,,j] ))
    for(j in 2:nt) {
      phi.ts[,j]<-rho.ts[j]*phi.ts[,j-1] + rmvnorm(1, rep(0, max.clust), lambda[j]*temp.H)
    }
    
    #Matrix of dynamic random effects
    for (j in 1:nt) {
      phi.effect[,j] <- phi.ts[,j][g.mat[,j]]
    }

    #Update lambda
    phi.temp <- cbind(phi.effect[,2 : (nt)], rep(0, nt))
    for (j in 1:max.clust) {
      lambda[j] <- rgamma(1, 1 + nt*ns/2, 1 + 0.5 * (t(rowSums(phi.effect - phi.temp*rho.ts[j])) %*% solve(H) %*% rowSums(phi.effect - phi.temp*rho.ts[j]) ) ) ^ (-1)
    }
  
    #update the rho.ts
    for (j in 1:max.clust) {
     v.temp <- 1/lambda[j] * t(phi.effect[,(nt-1)]) %*% solve(H) %*% phi.effect[,(nt-1)] + 10^(-4)
     d.temp <- 1/lambda[j] * t(phi.effect[,(nt-1)]) %*% solve(H) %*% phi.effect[,(nt-1)]
     rho.temp <- max(rnorm(1, v.temp^(-1) * d.temp, v.temp^(-1)),-1)
     rho.temp <- min(rho.temp,1)
     rho.ts[j] <- rho.temp
    }
    
    
    ##General updating/recording
        
    keeprho[i]<-rho
    keepbeta[i,]<-beta[,1]
    cluster.count[i] <- length(unique(g.mat[,1]))
    
    count<-count+1
    if(count==display){
      par(mfrow=c(1,1))
      plot(keeprho[1:i])
      plot(cluster.count[1:i])
      #sdmu<-(mu-min(mu))/(max(mu)-min(mu))
      plot(z[,1],z[,2],col=g.mat[,1], type = 'p', pch = 19)
      count<-0
    }
    
    if(i>burn){
      afterburn<-afterburn+1
      summu<-summu+mu.mat[g.mat[,1]]
      summu2<-summu2+mu.mat[g.mat[,1]]^2
      sumtp<-sumtp+probs[,n.terms]
    }
    
  }
  
  post.mn<-summu/afterburn
  post.sd<-summu2/afterburn-post.mn^2
  truncprob<-sumtp/afterburn
  
  list(truncprob=truncprob,beta=keepbeta[burn:runs,],
       rho=keeprho[burn:runs],post.mn=post.mn,post.sd=post.sd, membership = g.mat, mu = mu.mat[1,], knot.1 = knot.mat[,,1], dynamic.effect = phi.effect, dynamic.unique = phi.ts, probs.cluster1 = probs.cluster1)
}
