library(mvtnorm)
library(evd)

#Simulate data
nt<-50   #number of years
ns<-20   #number of spatial locations

#generate the spatial locations
S<-matrix(runif(2*ns),ns,2)
d<-as.matrix(dist(S,upper=T,diag=T))
covar<-exp(-d/.25)

#Generate the data
set.seed(0820)
y<-rmvnorm(nt,rep(0,ns),covar)
true_gev_location<-20*(S[,1]-.5)^2
true_gev_scale<-1
true_gev_shape<--.1
for(t in 1:nt){
  y[t,]<-qgev(pnorm(y[t,]),loc=true_gev_location,
              shape=true_gev_shape,
              scale=true_gev_scale) 
}
y<-t(y)
