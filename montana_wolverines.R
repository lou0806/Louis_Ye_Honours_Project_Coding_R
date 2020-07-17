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

gsdp_no_cluster(y,x,mx.siga=2, mx.sigb = 2,mx.taua = 2,mx.taub = 2)

SSB(y.ssb,z=x.ssb)
