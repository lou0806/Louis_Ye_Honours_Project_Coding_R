library(LatticeKrig) # quilt.plot
library(maps) # map()
library(Rcpp)
s_new <- read.csv("Data_rainfall/predlocs.dat")
set.seed(1)

load("Data_rainfall/Y.RData")

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
z <- bound.map(ylatlon)
y <- t(Y)

output <- SSB.mod(y =  matrix(y[,1],ncol = 1), z = z, x = x, n.terms = n, DOF=1,mx.sige=2,mx.sigs=2)
plot(z[,1],z[,2],col=output$membership, type = 'p', pch = 19)
points(clusters_dist, type = 'p')


## Functions

bound.map <- function(z) {
  minLong<-min(z[,1],na.rm = TRUE)
  maxLong<-max(z[,1],na.rm = TRUE)
  minLat<-min(z[,2],na.rm = TRUE)
  maxLat<-max(z[,2],na.rm = TRUE)
  
  longVec<-(z[,1] - minLong)/(maxLong - minLong)
  latVec<-(z[,2] - minLat)/(maxLat - minLat)
  
  return(cbind(longVec, latVec))
}

# lon,lat,val
viewPred <- function(x,latlon,main.plot='',bks=c(14,40)) {
  quilt.plot(latlon[,2],latlon[,1],x,
             fg='grey90',bty='n',main=main.plot,
             ylim=range(latlon[,1])+c(-1,1),
             xlim=range(latlon[,2])+c(-1,1),
             breaks=seq(bks[1],bks[2],len=length(x)+1),
             col= colorRampPalette(c('dark blue','grey90','dark red'))(length(x)))
  map('county',add=T,col='grey')
  map('state',add=T,col='grey60',lwd=2)
}


viewYearJuly <- function(yr,y,bks=seq(14,40,len=101)) {
  ind <- which(uyears==yr)
  quilt.plot(ylatlon[,2],ylatlon[,1],y[ind,],
             fg='grey90',bty='n',main=yr,
             ylim=range(ylatlon[,1])+c(-1,1),
             xlim=range(ylatlon[,2])+c(-1,1),
             breaks=bks,
             col= colorRampPalette(c('dark blue','grey90','dark red'))(100))
  map('county',add=T,col='grey')
  map('state',add=T,col='grey60',lwd=2)
}
viewYearJuly(1989,Y) #1985 - 2004
