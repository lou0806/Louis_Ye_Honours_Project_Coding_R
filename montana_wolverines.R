testdat2<- read.csv('Data/Montanawolverines.csv')

minLong<-min(testdat2$long,na.rm = TRUE)
maxLong<-max(testdat2$long,na.rm = TRUE)
minLat<-min(testdat2$lat,na.rm = TRUE)
maxLat<-max(testdat2$lat,na.rm = TRUE)

longVec<-(testdat2$long - minLong)/(maxLong - minLong)
latVec<-(testdat2$lat - minLat)/(maxLat - minLat)

modDat<-data.frame(testdat2[,1:10],longVec = longVec,latVec = latVec)

y<-as.matrix(testdat2[,1:10], ncol = 10)
x<-matrix(c(longVec,latVec),ncol = 2)

y.ssb <- matrix(unlist(testdat2[,1:10]))
x.ssb <- matrix(c(rep(x[,1],ncol(y)),rep(x[,2],ncol(y))), ncol = 2)

mx.allele<-6
y.spatClust<-array(0,c(nrow(x),ncol(y),mx.allele))
for(i in 1:nrow(x)){for(l in 1:ncol(y)){for(rep in 1:mx.allele){
  y.spatClust[i,l,rep]<-y[i,l]
}}}

gsdp_no_cluster(y,x,mx.siga=6, mx.sigb = 6,mx.taua = 6,mx.taub = 6)

SpatialClust(z=y.spatClust,n.all = rep(3,ncol(y)),x=x,spatial=T,genetic=T)

SSB(y.ssb,z=x.ssb)
