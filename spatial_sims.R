#Generated Data
library(ggplot2)

n=40 #No. of Locations
M=100 #No. of Observations
#A=2
n.clust=4 #"Correct" Number of Clusters
clusters_dist <- matrix(rep(c(.25,.5,.75),2),ncol = 2)
  #matrix(c(runif(n.clust,0,1),runif(n.clust,0,1)),ncol = 2)
clust_theta <- rnorm(n.clust,50,10)

#explore generated cluster effects
par(mfrow = c(1,1))
plot(clusters_dist)
print(clust_theta)

x<-matrix(runif(n*2,0,1),ncol = 2)
beta<-matrix(c(1,1),ncol = 2)
y<-matrix(rep(0,n*M),ncol = M)
tau<-.2

for (i in 1:n){
  temp.dex <- which.min((x[i,1] - clusters_dist[,1])^2 + (x[i,2] - clusters_dist[,2])^2)
  y[i,] = rnorm(M,mean=beta%*%x[i,] + clust_theta[temp.dex],sd=tau^2)
}

y.ssb <- matrix(unlist(y)[1:n])
x.ssb <- matrix(c(rep(x[,1],ncol(y)),rep(x[,2],ncol(y))), ncol = 2)[1:n,]

output <- SSB.mod(y = y.ssb, z = x.ssb)


##TODO: Vary alpha, not just the hyperparameters
y[1,]
loc.1 <- gsdp(y,x,mx.siga=1, mx.sigb = .25,mx.taua = 1,mx.taub = .25,a.alpha = 3,b.alpha = 1, loc = 1 , mx.bphi = .1)
  loc.2 <- gsdp(y,x,mx.siga=1, mx.sigb = .25,mx.taua = 1,mx.taub = .25,a.alpha = 10,b.alpha = 1, loc = 2, mx.bphi = .5)
loc.3 <- gsdp(y,x,mx.siga=1, mx.sigb = .25,mx.taua = 1,mx.taub = .25,a.alpha = 10,b.alpha = 1, loc = 3)


gsdp.cluster(y,x,mx.siga=1, mx.sigb = .25,mx.taua = 1,mx.taub = .25,a.alpha = 3,b.alpha = 1,  mx.bphi = .1)

points(clusters_dist,pch = 0)


SpatialClust(z=y.spatClust,n.all = rep(3,ncol(y)),x=x,spatial=T,genetic=T)
