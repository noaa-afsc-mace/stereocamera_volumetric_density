rm(list = ls())
library(plot3D)
library(ptinpoly)


load('data/pollock.rda')

# get the volume estimated
vol <- get_vol_func(pollock$cal,max_extent=7, grid_size=0.05, plotting=FALSE)

# define detection probability parameters
d50=3# this is the range at which 50% of criters are detection

sl=-1# this is the slope of detection curve, equivalent to selection range from trawl selectivity

scaling=1000 # appropriate for target positions in m
# loop starts here
n = 100000 # number of targets - in this case its targets over a set of frames, so much denser than what you would expect
reps=1000
estimated_density=numeric(length=reps)
for (i in 1:reps){

  # generate density data
  # now we add some random targets to the general space
  xf=runif(n)*10-5
  yf=runif(n)*10
  zf=runif(n)*10-5
  # find ones that are in both cones
  target_positions=matrix(c(xf,yf,zf),length(xf),3)
  targets_in_both=find_targets_in_view(pollock$cal,scaling, target_positions)

  # estimate probability of detection
  fish_range=sqrt(targets_in_both[,1]^2+targets_in_both[,2]^2+targets_in_both[,3]^2)
  prob=1./(1+9^((d50-fish_range)/sl))
  detected=numeric(length(prob))
  for (j in 1:length(prob)){detected[j]=rbinom(1,1,prob[j])}
  detected_fish_ranges=fish_range[which(detected==1)]
  detection_data=prep_detection_data(target_ranges=detected_fish_ranges,
                                     vol_func=vol$vol_func, nbins=NULL, method='max', nvals=2, loc_dens=NULL, plotting=FALSE)

  # round inputs to remove glm warning
  detection_data$exp_count <- round(detection_data$exp_count)
  detection_data$obs_count <- round(detection_data$obs_count)

  # fit detection function
  out <- fit_density_function(detection_data,method='logistic glm',formula=NULL, plotting=FALSE)


  # integrate
  eff_vol_pollock=integrate(eff_vol_func,lower =0,
                              upper =max(fish_range), vol$vol_func, out$detect.function)$value
  # estimate density
  estimated_density[i]=length(detected_fish_ranges)/eff_vol_pollock
}


known_density = n/(10^3)# this is in fish/m3
hist(estimated_density)
abline(v=known_density,col="red")
abline(v=mean(estimated_density),col="green")
# take home - we are "close-ish" with this method
# the change-in-volume function can be assumed to be "fixed" as it is dependent
# on physical properties of the camera system.  The detection function, and by
# extension the "effective volume", as it pertains to a specific species/imaging
# condition is what is needed to be estimated from existing data


