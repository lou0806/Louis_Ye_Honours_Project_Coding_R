#Adapted from https://www4.stat.ncsu.edu/~bjreich/code/SSB.R

library(MCMCpack)
library(distr)
library(mvtnorm)

library(LatticeKrig) # quilt.plot
library(maps) # map()
library(Rcpp)

#Compute v_k(s) (function from Reich(2007))
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

#take all the vs and compute probabilities (function from Reich(2007))
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

#generate categorical variables: (function from Reich(2007))
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

order.allocation <- function(mat) {
  nt <- ncol(mat)
  for (t in 1:nt) {
    temp <- 1
    for (i in sort(unique(mat[,t]))) {
      mat[,t][mat[,t] == i] <- temp
      temp <- temp + 1
    }
  }
  return(mat)
}

#The whole Shebang!
SSB.dynamic<-function(y,x=NA,z,DOF=1,mx.sige=1,mx.sigs=1,n.terms=100,
                      runs=1000,burn=300,display=10, testingRainfall = TRUE) {
  #y:       data
  #x:       covariates
  #z:       n x 2 matrix of coordinates, scaled to [0,1]
  #n.terms: number of terms in the mixture dist.
  #runs:    number of MCMC samples
  #burn:    number of samples to discard
  #display: how often to display results
  #DOF:     v~beta(1,DOF)
  #mx.sig:  the sds are U(0,mx.sd)
  #testingRainfall: Can remove this and anything associated, was quick fix for the density plots
  
  ns <- nrow(y) #number of locations
  nt <- ncol(y) #number of years/replications
  #n.terms <- ncol(y)
  
  n<-ns
  if(is.na(max(x))){x<-matrix(1,n,1)}
  p<-ncol(x)
  
  #Standardize the outcomes and predictors:
  if(min(z)<0 | max(z)>1){print("THE SPATIAL COORDINATES ARE NOT STANDARDIZED!")}   
  
  dens.1 <- dens.19 <- dens.34 <- rep(0, runs)
  
  #initial values
  beta.mat<-matrix(rep(0,nt*p), ncol = nt)
  cluster.count <- rep(0, runs)
  v <- rep(.9,n.terms) #originally .9, can change to influence initial clustering
  sige <- matrix(rep(mx.sige/2,n.terms*nt), nrow = nt)
  taue <- 1/sige^2
  mu.mat <- matrix(rep(0,n.terms*nt), nrow = nt)
  mu.value <- mu.mat
  sigs <- rep(mx.sigs/2, nt)
  taumu <- 1/sigs^2
  knot <- matrix(runif(2*n.terms,0,1),2,n.terms)
  rho <- rep(.5,nt) #influences clustering
  
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
      #COV[t]<-solve(t(x)%*%diag(taue[t,][g.mat[,t]])%*%x)
      #mn<-COV[t]%*%t(x)%*%diag(taue[g.mat[,t]])%*%(y[,t]-mu.mat[g.mat[,t]])
      beta.mat[,t] <- mean(y[,t])
      #beta.mat[,t] <- mn + t(chol(COV[t])) #ORIGINAL CODE FOR HIGHER DIMENSIONAL DATA (more covariates)
        #Altered for no covariate data
      r[,t]<-y[,t]-x%*%beta.mat[,t]
      knot <- knot.mat[,,t]
      old.g.mat <- g.mat
    
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
        MHrate<-dnorm(r[s,t],mu.mat[t,][cang[s]],sige[cang[s]],log=T)-
          dnorm(r[s,t],mu.mat[t,][g.mat[s,t]],sige[g.mat[s,t]],log=T)
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
        rho[t]<-canrho
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
            canvs[,k]<-onevs(z,rho[t],knot[,k],canv[k])
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
            canvs[,k] <- onevs(z,rho[t],canknot[,k],v[k])
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
      
      if(t == nt) {
        probs.last.cluster <- probs[,g.mat[1,]]
      }
    }
    
    #update Paci-Finazzi:
    #Keeping a track of knots (dynamic)
    
    #order g.mat for the Paci-Finazzi time series
    g.mat.paci <- order.allocation(g.mat)
    phi.ts.old <- phi.ts #Debug term
    n.clust <- rep(0, nt)
    max.clust <- max(max(length(unique(as.vector(g.mat.paci)))),2)
    max.clust.knot <- max(max(length(unique(as.vector(g.mat)))),2)
    
    knot.clust <- array(rep(knot, nt), dim = c(2, max.clust.knot, nt))
    phi.ts <- matrix(rep(0, max.clust.knot*nt), ncol = nt)
    phi.effect <- matrix(rep(0, nt*ns), ncol = nt)
    temp.length <- rep(0, nt)
    if(length(rho.ts) < max.clust) rho.ts <- c(rho.ts, rep(0, (max.clust - length(rho.ts))))
    
    for (t in 1:nt) {
      n.clust[t] <- length(unique(g.mat.paci[,t]))
      temp.length[t] <- max.clust.knot - n.clust[t]#ncol(knot.mat[,,t][,unique(g.mat[,t])])
      if (temp.length[t] > 0){
        temp.knot <- cbind(knot.mat[,,t][,unique(g.mat[,t])] 
                           , matrix(rep(0.5, (temp.length[t]*2)), nrow = 2))
      } else {
        temp.knot <- knot.mat[,,t][,unique(g.mat[,t])]
      }
      knot.clust[,,t] <- temp.knot
    }
    if (max(temp.length) > 0) {
      lambda <- c(lambda, rep(0.1, max(temp.length)))
    }

    #Unique random effects 1st time round
    phi.ts[,1] <- rmvnorm(1, rep(0, max.clust.knot), lambda[1]*exp.H(decay, t(knot.clust[,,1] )))
    for(k in 2:nt) {
      for (l in 1:max.clust){
        temp.H <- exp.H(decay, t(knot.clust[,,k] ))
        temp.mat <- rmvnorm(1, rep(0, max.clust.knot), lambda[l]*temp.H)
        phi.ts[l,k]<-rho.ts[k]*phi.ts[l,k-1] + temp.mat[l]
      }
    }
 
    #Matrix of dynamic random effects
    for (j in 1:nt) {
      phi.effect[,j] <- phi.ts[,j][g.mat.paci[,j]]
    }

    #Update lambda
    phi.temp <- cbind(phi.effect[,2 : (nt)], rep(0, nt))
    lambda.old <- lambda #Debug term
    lambda <- rep(0, max.clust)
    
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
    keeprho[i]<-rho[1]
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
      summu<-summu+mu.mat[g.mat[,nt]]
      summu2<-summu2+mu.mat[g.mat[,nt]]^2
      sumtp<-sumtp+probs[,n.terms]
    }
    
    print(i)
    
    #Density plots, can be deleted
    
    if (testingRainfall == TRUE) {
    
      mu.post.1 <- mu.mat[nt,][g.mat[1,nt]]
      mu.post.19 <- mu.mat[nt,][g.mat[19,nt]]
      mu.post.34<- mu.mat[1,][g.mat[34,1]]
      
      sd.pred.1 <- taue[nt,][g.mat[1,nt]]
      sd.pred.19 <- taue[nt,][g.mat[19,nt]]
      sd.pred.34 <- taue[1,][g.mat[34,1]]
      
      dens.1[i] <- rnorm(1, mean(y[,nt]) + mu.post.1, sd.pred.1)
      dens.19[i] <- rnorm(1, mean(y[,nt]) + mu.post.19, sd.pred.19)
      dens.34[i] <- rnorm(1, mean(y[,nt]) + mu.post.34, sd.pred.34)
      
      hist(dens.1, probability = TRUE)
      hist(dens.19, probability = TRUE)
      hist(dens.34, probability = TRUE)
    }
  }
  
  #Checking that the algorithm doesn't produce abnormal results
  post.mn<-summu/afterburn
  post.sd<-summu2/afterburn-post.mn^2
  truncprob<-sumtp/afterburn
  
  #Predictive surface
  surface.pred <- rep(0, nt)
  error.ts <- rmvnorm(1, rep(0, max.clust.knot), lambda[l]*temp.H)
  
  mn.post <- beta.mat[,nt]
  mu.post <- mu.mat[nt,][g.mat[,nt]]
  
  mu.pred <- mn.post*x + mu.post
  sd.pred <- taue[nt,][g.mat[,nt]]
  
  surface.pred <- mu.pred + phi.effect[,nt]*rho.ts[g.mat.paci[,nt]] + error.ts[g.mat.paci[,nt]]
  
  list(truncprob=truncprob,beta=keepbeta[burn:runs,],
       rho=keeprho[burn:runs],post.mn=post.mn,post.sd=post.sd, membership = g.mat, 
       mu = mu.mat[1,], knot.1 = knot.mat[,,1], dynamic.effect = phi.effect, 
       dynamic.unique = phi.ts, probs.cluster1 = probs.cluster1, rho.ts = rho.ts, 
       probs.last.cluster = probs.last.cluster, sd.pred = sd.pred, surface.pred = surface.pred, 
       ts.membership = g.mat.paci, mn.post = mn.post, mu.post = mu.post, 
       knot.pred = knot.mat[,,nt], lastprobs = probs, old.g.mat = old.g.mat,
       dens.1 = dens.1, dens.19 = dens.19, dens.34 = dens.34,
       error.ts = error.ts)
}

