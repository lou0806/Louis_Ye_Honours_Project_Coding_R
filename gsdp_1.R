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
  
  
}