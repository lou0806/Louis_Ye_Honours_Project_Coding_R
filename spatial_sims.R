#Generated Data

n=100 #No. of Locations
M=10 #No. of Observations
#A=2
n.clust=4 #"Correct" Number of Clusters
clusters_dist <- matrix(rep(c(.25,.5,.75),2),ncol = 2)
  #matrix(c(runif(n.clust,0,1),runif(n.clust,0,1)),ncol = 2)
clust_theta <- rnorm(n.clust,0,1)

#explore generated cluster effects
par(mfrow = c(1,1))
plot(clusters_dist)
print(clust_theta)

x<-matrix(runif(n*2,0,1),ncol = 2)
beta<-matrix(c(1,1),ncol = 2)
y<-matrix(rep(0,n*M),ncol = M)
tau<-.01

for (i in 1:n){
  temp.dex <- which.min((x[i,1] - clusters_dist[,1])^2 + (x[i,2] - clusters_dist[,2])^2)
  y[i,] = rnorm(M,mean=beta%*%x[i,] + clust_theta[temp.dex],sd=tau^2)
}

y.ssb <- matrix(unlist(y))
x.ssb <- matrix(c(rep(x[,1],ncol(y)),rep(x[,2],ncol(y))), ncol = 2)

##TODO: Vary alpha, not just the hyperparameters
gsdp(y,x,mx.siga=1, mx.sigb = .25,mx.taua = 1,mx.taub = .25,a.alpha = 1,b.alpha = 1)
points(clusters_dist,pch = 0)

SSB(y.ssb,z=x.ssb)

SpatialClust(z=y.spatClust,n.all = rep(3,ncol(y)),x=x,spatial=T,genetic=T)
