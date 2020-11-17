### SIMULATED DATA STUDY

library(ggplot2)
library(LatticeKrig) # quilt.plot
library(maps) # map()
library(Rcpp)
source("ssb_dynamic.R")

set.seed(2)

n = 100 #No. of Locations
M = 5 #No. of Observations
n.clust = 3 #"Correct" Number of Clusters

clusters_dist <- matrix(rep(seq(from = 0, to = 1, length.out = n.clust),2),ncol = 2)
clust_theta <- rnorm(n.clust,50,20)

#explore generated cluster effects
par(mfrow = c(1,1))
plot(clusters_dist)
print(clust_theta)

x<-matrix(runif(n*2,0,1),ncol = 2)
beta<-matrix(c(1,1),ncol = 2)
y<-matrix(rep(0,n*M),ncol = M)
tau<-1

time.effect <- seq(from = 0, to = 10, length.out = M)
for (i in 1:n){
  temp.dex <- which.min((x[i,1] - clusters_dist[,1])^2 + (x[i,2] - clusters_dist[,2])^2)
  y[i,] = rnorm(M,mean=beta%*%x[i,] + clust_theta[temp.dex],sd=tau^2) + time.effect
}

points(clusters_dist, type = 'p')

z <- x #coordinates
x <- matrix(rep(1, n), ncol = 1)

nt <- ncol(y)
ns <- nrow(y)
#g.true <- matrix(rep(0,ns*nt), ncol = nt)

## RUN THIS FUNCTION
output.dynamic <- SSB.dynamic(y =  y, z = z, x = x, n.terms = n, DOF=1, mx.sige=2, mx.sigs=2, runs = 100, burn = 50)

## Plots of Cluster Analysis
plot(z[,1],z[,2],col=output.dynamic$membership[,1], type = 'p', pch = 19)
quilt.plot(z[,1],z[,2],output.dynamic$probs.cluster1[,1],
           col= colorRampPalette(c('dark blue','grey90','dark red'))(100),
           xlim = c(0,1), ylim = c(0,1))
quilt.plot(z[,1],z[,2],output.dynamic$probs.last.cluster,
           col= colorRampPalette(c('dark blue','grey90','dark red'))(100),
           xlim = c(0,1), ylim = c(0,1))
plot(z[,1],z[,2],col=output.dynamic$membership[,5], type = 'p', pch = 19)

## Surface Prediction
surface.predicted <- output.dynamic$surface.pred
quilt.plot(z[,1],z[,2],surface.predicted,
           col= colorRampPalette(c('dark blue','grey90','dark red'))(100),
           xlim = c(0,1), ylim = c(0,1))

## Predictive Parameters for Smooth Surface Generation
rho.ts.pred <- output.dynamic$rho.ts[1:length(unique(output.dynamic$membership[,1]))]
effect.ts.pred <- output.dynamic$dynamic.effect[,M]
sd.pred <- output.dynamic$sd.pred

# Time series component
cluster.effect.pred <- output.dynamic$dynamic.unique
ts.membership.pred <- output.dynamic$ts.membership

# SSB Surface Component
mu.post.pred <- output.dynamic$mu.post
mn.post.pred <- output.dynamic$mn.post
sd.pred <- output.dynamic$sd.pred
membership.pred <- output.dynamic$membership

knots.pred <- output.dynamic$knot.pred

## Misc
output.dynamic$probs.cluster1
output.dynamic$dynamic.effect
output.dynamic$membership[,1]

output.dynamic$rho.ts[1:length(unique(output.dynamic$membership[,1]))]






