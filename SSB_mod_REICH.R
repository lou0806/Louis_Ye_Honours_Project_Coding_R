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


#The whole Shebang!
SSB.mod<-function(y,x=NA,z,DOF=1,mx.sige=1,mx.sigs=1,n.terms=100,
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
  beta <- rep(mean(y),p)
  v <- rep(.9,n.terms) #originally .9, changing to influence initial clustering
  sige <- rep(mx.sige/2,n.terms)
  taue <- 1/sige^2
  mu <- rep(0,n.terms)
  sigs <- mx.sigs/2
  taumu <- 1/sigs^2
  knot <- matrix(runif(2*n.terms,0,1),2,n.terms)
  rho <- .3 #originally .5, influences clustering by space
  
  g.mat <- matrix(rep(1,n*nt), ncol = nt) #membership
  g <- g.mat
  y.mat <- y
  
  vs<-matrix(0,n,n.terms)
  for(k in 1:n.terms){
    vs[,k]<-onevs(z,rho,knot[,k],v[k])
  }
  probs<-makeprobs(vs)
  
  sumtp<-summu<-summu2<-rep(0,n)
  count<-afterburn<-0
  keeprho<-rep(0,runs)
  keepbeta<-matrix(0,runs,p)
  keepmn <- rep(0, runs)
  keepmu <- rep(0, runs)
  keepr <- rep(0, runs)
  cluster.count <- rep(0, runs)
  
  for(i in 1:runs) {
    
    #Update beta
    COV<-solve(t(x)%*%diag(taue[g])%*%x)
    mn<-COV%*%t(x)%*%diag(taue[g])%*%(y-mu[g])
    beta <- mean(y) #use mean(y) when there are no covariates
      ##mean(y)
      ##mn + t(chol(COV))
    r<-y-x%*%beta
    
    for(j in 1:n.terms){
      #update mu
      nobs<-sum(g==j)
      mu[j]<-rnorm(1,0,sigs)
      if(nobs>0){
        mu[j]<-rnorm(1,mean(r[g==j]),1/sqrt(taumu+nobs*taue[j]))
      }
      
      #update sige
      cansige<-sige
      cansige[j]<-rnorm(1,sige[j],.1)
      if(sum(g==j)>0 & cansige[j]>0 & cansige[j]<mx.sige){
        MHrate<-sum(dnorm(r[g==j],mu[j],cansige[j],log=T)-
                      dnorm(r[g==j],mu[j],sige[j],log=T))
        if(runif(1,0,1)<exp(MHrate)){sige<-cansige} 
      }
      if(sum(g==j)==0 & cansige[j]>0 & cansige[j]<mx.sige){
        sige<-cansige 
      }
      taue<-1/sige^2
    }
    
    #update sigs:
    oldg <- g
    cansigs<-rnorm(1,sigs,.1)
    if(cansigs>0 & cansigs<mx.sigs){
      MHrate<-sum(dnorm(mu,0,cansigs,log=T)-
                    dnorm(mu,0,sigs,log=T))
      if(runif(1,0,1)<exp(MHrate)){sigs<-cansigs} 
    }
    taumu<-1/sigs^2
    
    #update g
    for(s in 1:n.terms){
      cang<-g
      cang[s]<-rcat(probs[s,])
      MHrate<-dnorm(r[s],mu[cang[s]],sige[cang[s]],log=T)-
        dnorm(r[s],mu[g[s]],sige[g[s]],log=T)
      if(runif(1,0,1)<exp(MHrate)){
        g<-cang
      } 
    }
    
    
    #update rho:
    canrho<-runif(1,0,1)
    canvs<-vs
    for(k in 1:n.terms){
      canvs[,k] <- onevs(z,canrho,knot[,k],v[k])
    }
    canprobs<-makeprobs(canvs)
    MHrate<-sum(log(canprobs[cbind(1:n,g)])-
                  log(probs[cbind(1:n,g)]))     
    if(runif(1,0,1)<exp(MHrate)){
      rho<-canrho;vs<-canvs;probs<-canprobs}
    
    #update v:
    for(k in 1:(n.terms-1)){
      if(max(g)<k){v[k]~rbeta(1,1,DOF)}
      if(max(g)>=k){
        canv<-v;canv[k]<-rnorm(1,v[k],.05)
        if(canv[k]>0 & canv[k]<1){
          canvs[,k]<-onevs(z,rho,knot[,k],canv[k])
          canprobs<-makeprobs(canvs)
          MHrate<-sum(log(canprobs[cbind(1:n,g)])-
                        log(probs[cbind(1:n,g)]))+
            log(dbeta(canv[k],1,DOF)/dbeta(v[k],1,DOF))
          if(runif(1,0,1)<exp(MHrate)){
            v<-canv;vs<-canvs
            probs<-canprobs
          }
        }
      }
    }
    
    #update knots:
    for(j in 1:2){
      for(k in 1:n.terms){
        canvs<-vs
        canknot<-knot
        canknot[j,k]<-rnorm(1,knot[j,k],.1)
        if(canknot[j,k]> 0 & canknot[j,k]<1){
          canvs[,k]<-onevs(z,rho,canknot[,k],v[k])
          canprobs<-makeprobs(canvs)
          MHrate<-sum(log(canprobs[cbind(1:n,g)])-
                        log(probs[cbind(1:n,g)]))     
          if(runif(1,0,1)<exp(MHrate)){
            knot<-canknot
            vs<-canvs;probs<-canprobs
          }
        }
      }
    } #Keep track of the log-likelihood
    
    keeprho[i]<-rho
    keepbeta[i,]<-beta
    keepmn[i] <- mn
    keepr[i] <- mean(r)
    cluster.count[i] <- length(unique(g))
    
    count<-count+1
    if(count==display){
      #par(mfrow=c(2,1))
      par(mfrow=c(1,1))
      #plot(keepbeta[,1])
      plot(keeprho[1:i])
      #plot(keepr[1:i])
      plot(cluster.count[1:i])
      sdmu<-(mu-min(mu))/(max(mu)-min(mu))
      plot(z[,1],z[,2],col=g, type = 'p', pch = 19)
      #plot(knot[,1],knot[,2])
      count<-0
    }
    
    if(i>burn){
      afterburn<-afterburn+1
      summu<-summu+mu[g]
      summu2<-summu2+mu[g]^2
      sumtp<-sumtp+probs[,n.terms]
    }
    
  }
  
  post.mn<-summu/afterburn
  post.sd<-summu2/afterburn-post.mn^2
  truncprob<-sumtp/afterburn
  
  list(truncprob=truncprob,beta=keepbeta[burn:runs,],
       rho=keeprho[burn:runs], residuals = keepr[burn:runs], post.mn=post.mn,post.sd=post.sd, 
       membership = g, knot = knot, mu = mu, taue = taue, sige = sige, keepmn = keepmn, 
       cluster.count = cluster.count, old.membership = oldg)
}

