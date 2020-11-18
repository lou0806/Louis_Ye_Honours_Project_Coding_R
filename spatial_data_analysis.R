## Data cleaning method adapted from http://luiarthur.github.io/sdp

library(LatticeKrig) # quilt.plot
library(maps) # map()
library(Rcpp)
s_new <- read.csv("Data_rainfall/predlocs.dat")
source("ssb_dynamic.R")

set.seed(1)

load("Data_rainfall/Y.RData")

## Bind the mapping to [0,1] x [0,1]

bound.map <- function(z) {
  minLong<-min(z[,1],na.rm = TRUE)
  maxLong<-max(z[,1],na.rm = TRUE)
  minLat<-min(z[,2],na.rm = TRUE)
  maxLat<-max(z[,2],na.rm = TRUE)
  
  longVec<-(z[,1] - minLong)/(maxLong - minLong)
  latVec<-(z[,2] - minLat)/(maxLat - minLat)
  
  return(cbind(longVec, latVec))
}

revert.map <- function(z, minLong, maxLong, minLat, maxLat) {

  longVec <- z[,1]*(maxLong - minLong) + minLong
  latVec <- z[,2]*(maxLat - minLat) + minLat
  
  return(cbind(longVec, latVec))
}

Ydim <- dim(Yout)
locs <- Ydim[1]
datcols <- Ydim[3]
ntimes <- Ydim[2]

uyears <- unique(Yout[1,,1])
TT <- length(uyears)
july <- which((1:ntimes)%%12%%7==0 & (1:ntimes)%%12!=0)
Y <- t(Yout[,july,5]) # Y is 20 x 100
# 20 years of july max temperatures
# 100 stations
#Y <- scale(Y,center=TRUE,scale=FALSE)[1:nrow(Y),]
ylatlon <- Yout[,1,3:4]
D <- as.matrix(dist(ylatlon))

# Y and ylatlon are the reponse/location covariates
minLong<-min(ylatlon[,1],na.rm = TRUE)
maxLong<-max(ylatlon[,1],na.rm = TRUE)
minLat<-min(ylatlon[,2],na.rm = TRUE)
maxLat<-max(ylatlon[,2],na.rm = TRUE)

z <- bound.map(ylatlon)
y <- t(Y)[,1:19]
pred.true <- t(Y)[,20]

output.rainfall <- SSB.dynamic(y =  y, z = z, x = x, n.terms = n, DOF=1,mx.sige=1,mx.sigs=1, runs = 100, burn = 50)

## Plots

plot(z[,1],z[,2],col=output.rainfall$membership[,1], type = 'p', pch = 19, main = "Clusters at t = 1")
plot(z[,1],z[,2],col=output.rainfall$membership[,5], type = 'p', pch = 19, main = "Clusters at t = 19")


first.clust.1 <- output.rainfall$membership[1,1]

viewPred.1 <- function() {
  quilt.plot(ylatlon[,2],ylatlon[,1],output.rainfall$probs.cluster1[,first.clust.1],fg='grey90',
           col= colorRampPalette(c('dark blue','grey90','dark red'))(100),
           nx=45,ny=45, main = "Pr of Belonging to Cluster 1 at t = 1")
  map('county',add=T,col='grey')
  map('state',add=T,col='grey60',lwd=2)
}
viewPred.1()

first.clust.19 <- output.rainfall$membership[19,1]
viewPred.last.cluster <- function() {
  quilt.plot(ylatlon[,2],ylatlon[,1],output.rainfall$probs.last.cluster[,first.clust.19],
           col= colorRampPalette(c('dark blue','grey90','dark red'))(100),
           nx=45,ny=45, main = "Pr of Belonging to Cluster 1 at t = 19")
  map('county',add=T,col='grey')
  map('state',add=T,col='grey60',lwd=2)
}

viewPred.last.cluster()



## Surface Prediction
surface.predicted <- output.rainfall$surface.pred
plot(density(output.rainfall$surface.pred), main = "Predicted Density vs True Density")
lines(density(pred.true), col = 'red')



viewPred.Surface <- function() {
  quilt.plot(ylatlon[,2],ylatlon[,1],surface.predicted,
           col= colorRampPalette(c('dark blue','grey90','dark red'))(100), main = "Predicted Surface",
           nx=45,ny=45)
  map('county',add=T,col='grey')
  map('state',add=T,col='grey60',lwd=2)
}
viewPred.Surface()
viewYearJuly(2004,Y) 


## Predictive Parameters for Smooth Surface Generation
rho.ts.pred <- output.rainfall$rho.ts[1:length(unique(output.rainfall$membership[,1]))]
effect.ts.pred <- output.rainfall$dynamic.effect[,M]
sd.pred <- output.rainfall$sd.pred

# Time series component
cluster.effect.pred <- output.rainfall$dynamic.unique
ts.membership.pred <- output.rainfall$ts.membership

# SSB Surface Component
mu.post.pred <- output.rainfall$mu.post
mn.post.pred <- output.rainfall$mn.post
theta.mean.pred <- mu.post.pred + mn.post.pred #MEANS
sd.pred <- output.rainfall$sd.pred
membership.pred <- output.rainfall$membership[,19]
hist(membership.pred, main = "Histogram of Cluster Membership", breaks = 20)

knots.pred <- output.rainfall$knot.pred
# Get knot locations
knots.pred.unique <- knots.pred[,unique(membership.pred)]
knots.unique.revert <- revert.map(z = t(knots.pred.unique), minLong = minLong, maxLong = maxLong, minLat = minLat, maxLat = maxLat)

# Visualisation of knot means and sd
viewPred.knotmeans <- function() {
  quilt.plot(knots.unique.revert[,2],knots.unique.revert[,1],unique(theta.mean.pred),
             col= colorRampPalette(c('dark blue','grey90','dark red'))(100), main = "Predicted Knot Means",
             nx=45,ny=45)
  map('county',add=T,col='grey')
  map('state',add=T,col='grey60',lwd=2)
}
viewPred.knotsd <- function() {
  quilt.plot(knots.unique.revert[,2],knots.unique.revert[,1],unique(sd.pred),
             col= colorRampPalette(c('dark blue','grey90','dark red'))(100), main = "Predicted Knot Sd",
             nx=45,ny=45)
  map('county',add=T,col='grey')
  map('state',add=T,col='grey60',lwd=2)
}
viewPred.knotmeans()
viewPred.knotsd()

# Density estimation
mu.effect <- unique(mu.post.pred)
which(unique(output.rainfall$membership[,19]) == temp.assign[1])

# Location 1
loc.1.density <- output.rainfall$dens.1
plot(density(loc.1.density[50:100]), main = "Density at Location 1")

#Location 19
loc.19.density <- output.rainfall$dens.19
plot(density(loc.19.density[50:100]), main = "Density at Location 19")
 

#Location 34
loc.34.density <- output.rainfall$dens.34
plot(density(loc.34.density[50:100]), main = "Density at Location 34")

#unique(output.rainfall$membership[,19])
#unique(mu.post.pred)










## Visualisation Functions (from Lui, Arthur)

# lon,lat,val
viewPred <- function(x,latlon,main.plot='',bks=c(14,40)) {
  quilt.plot(latlon[,2],latlon[,1],x,
             fg='grey90',bty='n',main=main.plot,
             ylim=range(latlon[,1])+c(-1,1),
             xlim=range(latlon[,2])+c(-1,1),
             breaks=seq(bks[1],bks[2],len=length(x)+1),
             col= colorRampPalette(c('dark blue','grey90','dark red'))(length(x)),
             nx=45,ny=45)
  map('county',add=T,col='grey')
  map('state',add=T,col='grey60',lwd=2)
}


viewYearJuly <- function(yr,y,bks=seq(14,40,len=101)) {
  ind <- which(uyears==yr)
  quilt.plot(ylatlon[,2],ylatlon[,1],y[ind,],
             fg='grey90',main=yr,
             col= colorRampPalette(c('dark blue','grey90','dark red'))(100),
             nx=45,ny=45)
  map('county',add=T,col='grey')
  map('state',add=T,col='grey60',lwd=2)
}
viewYearJuly(1985,Y) #1985 - 2004
viewYearJuly(2004,Y) 
