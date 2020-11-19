##https://www4.stat.ncsu.edu/~bjreich/software Reich BJ, Bondell HD (2011). A spatial Dirichlet process mixture model for clustering population genetic data. Biometrics


library(MCMCpack)

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

#update the stick-breaking parameters, v
new.v<-function(g,b,m){
   v<-rep(1,m)
   for(k in 1:(m-1)){
     v[k]<-rbeta(1,sum(g==k)+1,sum(g>k)+b)
   }
v}


#update the DP hyperparameter, b
new.b<-function(v,a,b){
   v[v>0.9999]<-0.9999
rgamma(1,sum(v>0)+a,b-sum(log(1-v)))}


#Generage a truncate normal   
rtnorm<-function(n,mu,sigma,L,U){
   l.prob<-pnorm(L,mu,sigma)
   u.prob<-pnorm(U,mu,sigma)
qnorm(runif(n,l.prob,u.prob),mu,sigma)}


#Generate a truncated gamma:
rtgamma<-function(n,a,b,L,U){
   l.prob<-pgamma(L,a,b)
   u.prob<-pgamma(U,a,b)
qgamma(runif(n,l.prob,u.prob),a,b)}


#The prior for rho.star
prior.rho.star<-function(rho,p0=0.5,lambda=c(0.001,20)){
    log((1-p0)*dexp(rho,1/lambda[1]) + p0*dunif(rho,0,lambda[2]))
}

rho.star.2.rho<-function(rho.star){1/rho.star}

# The log dirichlet density:
logdirichlet<-function(x,a,fudge=exp(-10)){
   x[x<fudge]<-fudge
   a[a<fudge]<-fudge
sum((a-1)*log(x)) + lgamma(sum(a))-sum(lgamma(a))}


#THIS IS THE MAIN MCMC FUNCTION!
SpatialClust<-function(z,n.all,x,
              m1=25,m2=25,
              spatial=T,genetic=T,
              runs=10000,burn=1000,update=10,
              adj.lam=T,
              common.range=T,
              p0=0.5,max.lambda=0.0001,n.rep=2,
              b.a=1,b.b=4,
              h.a=1,h.b=1,
              tau.a=.1,tau.b=.1){

    library(fields)


    #Determine the dimensions:
    n<-dim(z)[1]
    nloci<-dim(z)[2]
    N.all<-dim(z)[3]


    #Set the initial values:
    zbar<-matrix(0,N.all,nloci)
    missing<-is.na(z[,,1])    
    for(l in 1:nloci){
     for(j in 1:n.all[l]){
      zbar[j,l]<-mean(z[!missing[,l],l,j])/n.rep
     }
     zbar[1:n.all[l],l]<-0.001+zbar[1:n.all[l],l]
     zbar[1:n.all[l],l]<-zbar[1:n.all[l],l]/sum(zbar[1:n.all[l],l])
    }

    theta<-array(0,c(N.all,nloci,m1))
    for(l in 1:nloci){theta[1:n.all[l],l,]<-zbar[1:n.all[l],l]}
    v<-c(rep(0.5,m1-1),1)
    probs.g<-makeprobs(v)
    g<-as.vector(rcat(n,probs.g))

    rho.star<-rep(.2,nloci)
    mn.theta<-zbar*rho.star.2.rho(rho.star[1])
    
    mu1<-matrix(runif(m1*m2,-1,1),m1,m2)
    mu2<-matrix(runif(m1*m2,-1,1),m1,m2)
    u<-cbind(matrix(.5,m1,m2-1),1)
    probs.h<-u
    for(j in 1:m1){probs.h[j,]<-makeprobs(u[j,])}   
    h<-rep(1,n)
    b1<-b2<-1
    tau<-50+0*mu1;sig<-1/sqrt(tau)
    lambda<-cbind(0.001+0*n.all,20+0*n.all)
    if(adj.lam){
      lambda<-cbind(max.lambda*(n.all-1),20+0*5*(n.all-1))
    }

    #Keep track of stuff:

    keep.g<-matrix(0,runs,n)
    keep.mn.theta<-array(0,c(runs,N.all,nloci))
    keep.b1<-keep.b2<-keep.sig<-clusters<-clusters2<-keep.lambda<-rep(0,runs)
    keeprho.star<-matrix(0,runs,nloci)


    #Start the MCMC!
    for(i in 1:runs){

       #Fill in missings:
       if(sum(missing)>0){for(s in 1:n){for(l in 1:nloci){
          if(missing[s,l]){
             z[s,l,1:n.all[l]]<-
              rmultinom(1,n.rep,theta[1:n.all[l],l,g[s]])
          }
       }}}

       ###LIKELIHOOD:
       if(genetic){for(l in 1:nloci){for(j in 1:m1){
          zzz<-mn.theta[,l]
          if(sum(g==j)==1){zzz<-zzz+z[g==j,l,]}
          if(sum(g==j)>1){zzz<-zzz+apply(z[g==j,l,],2,sum)}
          theta[1:n.all[l],l,j]<-rdirichlet(1,zzz[1:n.all[l]])
       }}}

       if(genetic){
         for(l in 1:nloci){
          fudge<-50
          canrho.star<-rho.star
          canrho.star[l]<-rbeta(1,rho.star[l]*fudge,(1-rho.star[l])*fudge)
          if(canrho.star[l]>0.0000001){
           canmn.theta<-zbar[,l]*rho.star.2.rho(canrho.star[l])
           R<-prior.rho.star(canrho.star[l],p0=p0,lambda=lambda[l,])-
              prior.rho.star(rho.star[l],p0=p0,lambda=lambda[l,])+
              dbeta(rho.star[l],fudge*canrho.star[l],fudge*(1-canrho.star[l]),log=T)-
              dbeta(canrho.star[l],fudge*rho.star[l],fudge*(1-rho.star[l]),log=T)
           for(k in 1:m1){
           R<-R+logdirichlet(theta[1:n.all[l],l,k],canmn.theta[1:n.all[l]])-
                logdirichlet(theta[1:n.all[l],l,k],mn.theta[1:n.all[l],l])
           }       
           if(!is.na(exp(R))){if(runif(1)<exp(R)){
              mn.theta[,l]<-canmn.theta;rho.star<-canrho.star}}
          }
         }
      }

       ### DP:
       v<-new.v(g,b1,m1)
       probs.g<-makeprobs(v)
       b1<-new.b(v[-m1],b.a,b.b)
   
       #new g/h
       if(!spatial & !genetic){g<-rcat(n,probs.g)}
       if(!spatial & genetic){for(j in 1:n){
         ppp<-log(probs.g)
         for(l in 1:nloci){
            ppp<-ppp+as.vector(z[j,l,1:n.all[l]]%*%log(theta[1:n.all[l],l,]))
         }
         ppp<-exp(ppp-max(ppp))
         if(!is.na(sum(ppp))){if(min(ppp)>0){g[j]<-rcat(1,ppp)}}
       }}
       if(spatial){for(j in 1:n){
         ppp<-log(probs.h)+matrix(log(probs.g),m1,m2,byrow=F)-
              0.5*tau*(x[j,1]-mu1)^2-0.5*tau*(x[j,2]-mu2)^2
         if(genetic){for(l in 1:nloci){
             junk<-as.vector(z[j,l,1:n.all[l]]%*%log(theta[1:n.all[l],l,]))
             ppp<-ppp+matrix(junk,m1,m2,byrow=F)
         }}
         ppp<-exp(ppp-max(ppp))
         if(!is.na(sum(ppp))){if(min(ppp)>0){
           g[j]<-rcat(1,apply(ppp,1,sum))
           h[j]<-rcat(1,ppp[g[j],])
         }}
       }}

       if(spatial){for(j in 1:m1){for(k in 1:m2){
          xxx1<-x[g==j & h==k,1]
          xxx2<-x[g==j & h==k,2]
          VAR<-1/(.001+tau[j,k]*length(xxx1))
          mu1[j,k]<-rtnorm(1,tau[j,k]*VAR*sum(xxx1),sqrt(VAR),-1,1)
          mu2[j,k]<-rtnorm(1,tau[j,k]*VAR*sum(xxx2),sqrt(VAR),-1,1)
       }}}

       if(spatial & common.range){
         S1<-(x[,1]-mu1[cbind(g,h)])^2
         S2<-(x[,2]-mu2[cbind(g,h)])^2
         tau<-tau*0+rgamma(1,n+tau.a,sum(S1+S2)/2+tau.b)
         sig<-1/sqrt(tau)
       }
       if(spatial & !common.range){
         SSS<-(x[,1]-mu1[cbind(g,h)])^2+(x[,2]-mu2[cbind(g,h)])^2
         for(j in 1:m1){
           tau[j,]<-rgamma(1,sum(g==j)+tau.a,sum(SSS[g==j])/2+tau.b)
         }
         sig<-1/sqrt(tau)
       }


      #new u
      if(spatial){
         for(j in 1:m1){
           u[j,]<-new.v(h[g==j],b2,m2)
           probs.h[j,]<-makeprobs(u[j,])
         }
         b2<-new.b(u[,-m2],h.a,h.b)
      }

      keep.g[i,]<-g
      keep.b1[i]<-b1
      keep.b2[i]<-b2
      keep.sig[i]<-sig[1,1]
      keeprho.star[i,]<-rho.star
      keep.lambda[i]<-lambda[1,1]
      clusters2[i]<-sum(table(g)>1)
      clusters[i]<-sum(table(g)>0)

      if(i%%update==0){
        par(mfrow=c(1,2))
        plot(x,col=g,pch=19,main=paste("Cluster labels, iteration",i))
        plot(clusters2[1:i],type="l")
      }
    }

    pequal<-matrix(0,n,n)
    for(w1 in 1:n){for(w2 in 1:n){
      pequal[w1,w2]<-mean(keep.g[burn:runs,w1]==keep.g[burn:runs,w2])
    }}

    membership<-matrix(NA,n,m1)
    clusters<-clusters[burn:runs]
    clusters2<-clusters2[burn:runs]
    
    for(j in 1:m1){if(sum(clusters2==j)>10){  
      ggg<-keep.g[clusters2==j,]
      equals<-diag(n)
      for(j1 in 1:n){for(j2 in 1:n){
        equals[j1,j2]<-mean(ggg[,j1]==ggg[,j2])
      }}
      membership[,j]<-cutree(hclust(as.dist(1-equals)),k=j)
    }}

list(g=keep.g[burn:runs,],
     rho.star=keeprho.star[burn:runs,],
     clusters=clusters,
     clusters2=clusters2,
     pequal=pequal,
     membership=membership
)}



