



rm(list = ls())

library(likelihoodExplore)
library(tidyverse)
library(optimx)
targets<-read.csv("targets.csv",header=TRUE)
targets=targets %>% filter(SPECIES_GROUP=="herring")
targets$RANGE=targets$RANGE/100# change units from cm to m

load("pelagicam_2023_vol.RData")# this is our specific camera volume data, including the volume function

# here we bin the data for fitting the model.  This involves estimation the volume 
# for each range bin as well as how many animals are in the bin.
nbins=15
bin_edges=seq(min(targets$RANGE),max(targets$RANGE),length.out=nbins+1)
bin_mids=bin_edges[1:nbins]+(bin_edges[2]-bin_edges[1])/2

vol=numeric(length=length(bin_mids))
bin_counts=numeric(length=length(bin_mids))
for (i in 1:length(bin_mids)){
  v=integrate(vol_func,p_vol,lower = bin_edges[i], upper = bin_edges[i+1])
  vol[i]=v$value
  ind=which(targets$RANGE>bin_edges[i] & targets$RANGE<=bin_edges[i+1])
  bin_counts[i]=length(ind)
}

dens=bin_counts/vol
plot(bin_mids,dens,xlab="range (m)",ylab="density (fish/m^3)")

# logistic
detect_single_logistic = function(parm,vol,rng,fishcount){
  like=numeric(length=length(rng))
  for (i in 1:length(rng)){
    # expected number of fish per m3
    d_hat=parm[1]/(1+9^((parm[2]-rng[i])/parm[3]))
    # poisson likelihood?
    #ll=likpois(round(fishcount[i]/vol[i]), round(d_hat), log = TRUE)*fishcount[i]*-1
    # normal?
    ll=liknorm(fishcount[i]/vol[i], d_hat, log = TRUE)*fishcount[i]*-1 
    # key part here is to multiply the likelihood by the number of observations in the bin.
    # this way bins with more fish in them drive the fit
    like[i]=ll
  }
  return(sum(like))
}
# minimization
out=opm(par=c(25,4,-2), fn=detect_single_logistic, gr = NULL, vol=vol,
        rng=bin_mids,fishcount=bin_counts,method="Nelder-Mead")

disp_rng=seq(min(bin_mids),max(bin_mids),length.out = 100)
d_hat=out$p1/(1+9^((out$p2-disp_rng)/out$p3))
lines(disp_rng,d_hat,col="red")



# double logistic
detect_double_logistic = function(parm,vol,rng,fishcount){
  like=numeric(length=length(rng))
  for (i in 1:length(rng)){
    # expected number of fish per m3
    # this functin has a scale (param 1) and an rising and descending logistic component
    d_hat=(parm[1]*(1/(1+9^((parm[2]-rng[i])/parm[3])))*(1/(1+9^((parm[4]-rng[i])/parm[5]))))
    # poisson likelihood?
    #ll=likpois(round(fishcount[i]/vol[i]), round(d_hat), log = TRUE)*fishcount[i]*-1
    # normal?
    ll=liknorm(fishcount[i]/vol[i], d_hat, log = TRUE)*fishcount[i]*-1 
    # key part here is to multiply the likelihood by the number of observations in the bin.
    # this way bins with more fish in them drive the fit
    like[i]=ll
  }
  return(sum(like))
}


# minimization
out2=opm(par=c(25,1.8,2,3,-1), fn=detect_double_logistic, gr = NULL, vol=vol,
         rng=bin_mids,fishcount=bin_counts,method="Nelder-Mead")
#plot
d_hat=out2$p1*(1/(1+9^((out2$p2-disp_rng)/out2$p3)))*(1/(1+9^((out2$p4-disp_rng)/out2$p5)))
lines(disp_rng,d_hat,col="green")
