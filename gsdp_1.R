gsdp<-function(y,x,
                 spatial=T,varyparam=rep(F,3),nk=5,
                 pri.mn=rep(0,3),pri.sd=rep(100,3),
                 as=rep(0.1,3),bs=rep(0.1,3),mx.rho=10,
                 iters=5000,burn=1000,verbose=10,thin=1){
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
  #pri.mn,    the prior mean and standard deviation of the
  #pri.sd:        spatial mean GEV parameters
  ns<- nrow(y) #number of locations
  nt<- ncol(y) #number of years/replications
  
  #initial values
  beta<-rep(0,p)
  v<-rep(.9,n.terms)
  sige<-rep(mx.sige/2,n.terms);taue<-1/sige^2
  mu<-rep(0,n.terms)
  sigs<-mx.sigs/2;taumu<-1/sigs^2
  knot<-matrix(runif(2*n.terms,0,1),2,n.terms)
  rho<-.5
  g<-rep(1,n)
  vs<-matrix(0,n,n.terms)
  for(k in 1:n.terms){vs[,k]<-onevs(z,rho,knot[,k],v[k])}
  probs<-makeprobs(vs)
  
  v<-rep(.5,nk);v[nk]<-1
  D<-1
  probs<-makeprobs(v)
  g<-rep(1,nt)
  MU<-rep(0,nk)
  MU.prec<-1
  rho<-rep(mean(d),nk)
  
}