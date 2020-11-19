##https://www4.stat.ncsu.edu/~bjreich/code/SpatClustHelp.pdf

n=100
M=10
A=2
n.clust=3

#Allele frequencies for each cluster:
all.freq=array(1/A,c(M,A,n.clust))
all.freq[1,,1]=all.freq[2,,1]=all.freq[3,,1]=all.freq[4,,1]=c(.2,.8)
all.freq[1,,2]=all.freq[2,,2]=all.freq[3,,2]=all.freq[4,,2]=c(.8,.2)

#Spatial distribution of each cluster
s1=s2=c(-.5,.5,0)
bw=0.25

#Generate the genetic data:
n.all=rep(A,M)
z=array(0,c(n,M,A))
g=sample(1:n.clust,n,replace=T)
x=cbind(rnorm(n,s1[g],bw), rnorm(n,s2[g],bw))
for(i in 1:n){for(l in 1:M){for(rep in 1:2){
  samp<-sample(1:A,1, prob=all.freq[l,,g[i]])
  z[i,l,samp]= z[i,l,samp]+1
}}}

#Calculate the prior on the number of clusters:
fit=SpatialClust(z=z,n.all=n.all,x=x,spatial=F,genetic=F)
table(fit$clusters2)/length(fit$clusters2)

#Calculate the posterior of the number of clusters:
fit=SpatialClust(z=z,n.all=n.all,x=x)
table(fit$clusters2)/length(fit$clusters2)

#Plot the membership classification assuming three clusters:
plot(x[,1],x[,2],pch=19,col=fit$membership[,3])
