
output <- SSB.mod(y =  matrix(y[,1],ncol = 1), z = z, x = x, n.terms = n, DOF=1,mx.sige=1,mx.sigs=1)

points(clusters_dist, type = 'p')

plot(z[,1],z[,2],col=output$membership, type = 'p', pch = 19)
points(clusters_dist, type = 'p')
#unique(output$membership)
#points(t(output$knot[,unique(output$membership)]), type = 'p', pch =5)

output <- SSB.mod(y =  matrix(y[,3],ncol = 1), z = z, x = x, n.terms = 20, DOF=1,mx.sige=2,mx.sigs=2)

##Notes 15/10/2020:
## n.terms heavily influences clustering probabilities (away from first cluster)
# Idea: record the changes in allocation

##Notes 16/10/2020:
## Look at predictive differences with the Gelfand code (gsdp)

#Clear 2-cluster setup

#z <- rbind(cbind(rnorm(100,0,0.5), rnorm(100,0,0.5)), cbind(rnorm(150,5,0.5), rnorm(150,5,0.5)))
#z <- rbind(z, cbind(runif(100,min(z),max(z)), runif(100,min(z),max(z))))

for(t in 1:nt){
  y.vec <- matrix(y[,t],ncol = 1)
  output <- SSB.mod(y = y.vec, z = z)
  g.mat[,t] <- output$membership
  print(t)
}

for (t in 1:nt) {
  temp <- 1
  for (i in unique(g.mat[,t])) {
    g.mat[,t][g.mat[,t] == i] <- temp
    temp <- temp + 1
  }
  g.mat[,t]
}
#y.ssb <- matrix(unlist(y))
#x.ssb <- matrix(c(rep(x[,1],ncol(y)),rep(x[,2],ncol(y))), ncol = 2)
#
#output <- SSB.mod(y = y.ssb, z = x.ssb)







##TODO: Vary alpha, not just the hyperparameters
y[1,]
loc.1 <- gsdp(y,z,mx.siga=1, mx.sigb = .25,mx.taua = 1,mx.taub = .25,a.alpha = 3,b.alpha = 1, loc = 1 , mx.bphi = .1)
loc.2 <- gsdp(y,x,mx.siga=1, mx.sigb = .25,mx.taua = 1,mx.taub = .25,a.alpha = 10,b.alpha = 1, loc = 2, mx.bphi = .5)
loc.3 <- gsdp(y,x,mx.siga=1, mx.sigb = .25,mx.taua = 1,mx.taub = .25,a.alpha = 10,b.alpha = 1, loc = 3)


gsdp.cluster(y,x,mx.siga=1, mx.sigb = .25,mx.taua = 1,mx.taub = .25,a.alpha = 3,b.alpha = 1,  mx.bphi = .1)


SpatialClust(z=y.spatClust,n.all = rep(3,ncol(y)),x=x,spatial=T,genetic=T)