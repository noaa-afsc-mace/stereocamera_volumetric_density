rm(list = ls())
library(plot3D)
library(ptinpoly)


# now to run the simulation
n = 1000 # number of fish
reps=100
bins=vector(mode="numeric",length=10^3)
outmat=matrix(data=NA,nrow=reps,ncol=10^3)
plot(1, type = "n", xlim = c(0, 1000), ylim = c(0, 10),main = 'count per m3')
for (i in 1:reps){
  # now we add some random fish
  xf=runif(n)*10
  yf=runif(n)*10
  zf=runif(n)*10
  counter=1
  for(j in 1:10){
    for (k in 1:10){
      for (l in 1:10){

        bins[counter]=length(which(xf>j-1 & xf<=j & yf>k-1 & yf<=k & zf>l-1 & zf<=l))
        counter=counter+1
      }
    }
  }
  outmat[i,]=bins
  lines(bins, type = "l", col = "red")
}
