
rm(list = ls())
library(likelihoodExplore)
library(tidyverse)
library(optimx)
targets<-read.csv("detection_function/targets.csv",header=TRUE)
targets$RANGE=targets$RANGE/100# change units from cm to m
targets$HX=targets$HX/100
targets$HY=targets$HY/100
targets$HZ=targets$HZ/100

# for this version, we're implementing Cole's thought experiment where we constrain 
# the targets to a single non-expanding volume (like a fixed witch/height polygon 
# going away fropm the camera).  We lose a lot of data but do not have to worry 
# about the expanding volume


# here we bin the data for fitting the model. 
 
# these are the limits in m of the fish we are looking at +/- from 0(central axis)
x_limit=0.3
z_limit=0.3

# filter only to targets we want
targets=targets %>% filter(SPECIES_GROUP=="Adult_pollock" & HX>=-x_limit & HX<=x_limit & 
                             HZ>=-z_limit & HZ<=z_limit)


nbins=10
bin_edges=seq(min(targets$RANGE),max(targets$RANGE),length.out=nbins+1)
bin_mids=bin_edges[1:nbins]+(bin_edges[2]-bin_edges[1])/2
bin_counts=numeric(length=length(bin_mids))
for (i in 1:length(bin_mids)){
  ind=which(targets$RANGE>bin_edges[i] & targets$RANGE<=bin_edges[i+1])
  bin_counts[i]=length(ind)
}
# volume for each bin is "fixed", e.g. not range dependent
block_vol=2*x_limit*2*z_limit*(bin_edges[2]-bin_edges[1])

dens=bin_counts/block_vol
plot(bin_mids,dens,xlab="range (m)",ylab="density (fish/m^3)")

# logistic
detect_single_logistic = function(parm,rng,dens,fishcount){
  like=numeric(length=length(rng))
  for (i in 1:length(rng)){
    # expected number of fish per m3
    d_hat=parm[1]/(1+9^((parm[2]-rng[i])/parm[3]))
    # poisson likelihood?
    # ll=likpois(round(fishcount[i]/vol[i]), round(d_hat), log = TRUE)*fishcount[i]*-1
    # normal?
    ll=liknorm(dens[i], d_hat, log = TRUE)*fishcount[i]*-1 
    # key part here is to multiply the likelihood by the number of observations in the bin.
    # this way bins with more fish in them drive the fit
    like[i]=ll
  }
  return(sum(like))
}
# minimization
out=opm(par=c(500,3,-2), fn=detect_single_logistic, gr = NULL, dens=dens,
        rng=bin_mids,fishcount=bin_counts,method="Nelder-Mead")
#plot
disp_rng=seq(min(bin_mids),max(bin_mids),length.out = 100)
d_hat=out$p1/(1+9^((out$p2-disp_rng)/out$p3))
lines(disp_rng,d_hat,col="red")


# 
# # double logistic
# detect_double_logistic = function(parm,vol,rng,fishcount){
#   like=numeric(length=length(rng))
#   for (i in 1:length(rng)){
#     # expected number of fish per m3
#     # this functin has a scale (param 1) and an rising and descending logistic component
#     d_hat=(parm[1]*(1/(1+9^((parm[2]-rng[i])/parm[3])))*(1/(1+9^((parm[4]-rng[i])/parm[5]))))
#     # poisson likelihood?
#     #ll=likpois(round(fishcount[i]/vol[i]), round(d_hat), log = TRUE)*fishcount[i]*-1
#     # normal?
#     ll=liknorm(fishcount[i]/vol[i], d_hat, log = TRUE)*fishcount[i]*-1 
#     # key part here is to multiply the likelihood by the number of observations in the bin.
#     # this way bins with more fish in them drive the fit
#     like[i]=ll
#   }
#   return(sum(like))
# }
# 
# 
# # minimization
# out2=opm(par=c(700,1.5,0.5,3,-1), fn=detect_double_logistic, gr = NULL, vol=vol,
#          rng=bin_mids,fishcount=bin_counts,method="Nelder-Mead")
# #plot
# d_hat=out2$p1*(1/(1+9^((out2$p2-disp_rng)/out2$p3)))*(1/(1+9^((out2$p4-disp_rng)/out2$p5)))
# lines(disp_rng,d_hat,col="green")
# legend(x="topright",lty=c(1,1,1),col=c("black","red","green"), legend=c("data","single logistic","double logistic"))