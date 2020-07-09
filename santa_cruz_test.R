testdat<-X646985

#normalize data to [0,1]
#minLong<-min(testdat$LONGITUDE,na.rm = TRUE)
#maxLong<-max(testdat$LONGITUDE,na.rm = TRUE)
#minLat<-min(testdat$LATITUDE,na.rm = TRUE)
#maxLat<-max(testdat$LATITUDE,na.rm = TRUE)

#longVec<-(testdat$LONGITUDE - minLong)/(maxLong - minLong)
#latVec<-(testdat$LATITUDE - minLat)/(maxLat - minLat)


modDat<-data.frame(station = testdat$STATION,date = testdat$DATE,longVec = testdat$LONGITUDE,latVec = testdat$LATITUDE,TPCP = testdat$TPCP)

nRep<-length(unique(modDat$date))
nLoc<-length(unique(testdat$STATION))
stations<-unique(modDat$station)
dates<-unique(modDat$date)
responseMat<-data.frame(stations,matrix(rep(rep(NA,nLoc),nRep),ncol=nRep))

##Generates Y
##TODO: OPTIMISE
for (i in 1:nRep) {
  tempDat<-data.frame(station = modDat$station[dates[i]==modDat$date],TPCP = modDat$TPCP[dates[i]==modDat$date])
  for (j in 1:nLoc) {
    responseMat[j,i+1]<-tempDat$TPCP[tempDat$station==responseMat[j,1]][1]
  }
}
y<-matrix(as.numeric(unlist(responseMat[,1:nRep+1])),ncol=nRep,byrow = FALSE)

##Generates covariates X
##TODO: don't use mean to account for NA terms, fix later
x<-matrix(rep(NA,2*nLoc),ncol=2)
for (i in 1:nLoc){
  x[i,1]<-mean(modDat$longVec[stations[i]==modDat$station],na.rm=TRUE)
  x[i,2]<-mean(modDat$latVec[stations[i]==modDat$station],na.rm=TRUE)
}

#gsdp_no_cluster<-function(y,x,spatial=T,varyparam=rep(F,3),nk=5,iters=1000,burn=100,verbose=10,thin=1,mx.siga,mx.sigb,mx.taua,mx.taub, noClusters = NA)
gsdp_no_cluster(y,x,mx.siga=2, mx.sigb = 2,mx.taua = 2,mx.taub = 2)
