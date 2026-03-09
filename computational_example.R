# this is an example script with all the components in a single run, with the goal of
# providing a baseline for comparing a more structured approach - without plotting!

rm(list = ls())
library(StereoCamVolume)
runjellyfish=FALSE
runpollock=FALSE
runkrill=FALSE
runyelloweye=TRUE

if (runjellyfish==TRUE){
##########################################################################
################## Jellyfish (Sarsia sp.) example - Bering Sea 2024 ######

# this data set comes with the package now - replaced with a read function bc I don't know how this runs :-(
load('data/jellyfish.rda')

# get the volume estimated
vol <- get_vol_func(jellyfish$cal,max_extent=5, grid_size=0.05, plotting=TRUE, units='m3', seafloor_position=c(30,15,1))

# get density data
detection_data=prep_detection_data(target_ranges=jellyfish$targets$RANGE,
    vol_func=vol$vol_func, nbins=NULL, method='mean', nvals=3, loc_dens=NULL, plotting=TRUE)

# round inputs to remove glm warning
detection_data$exp_count <- round(detection_data$exp_count)
detection_data$obs_count <- round(detection_data$obs_count)

# fit detection function
out <- fit_density_function(detection_data,method='logistic glm',formula=NULL, plotting=TRUE)


# integrate
eff_vol_jellyfish=integrate(eff_vol_func,lower =0,
                            upper =max(jellyfish$targets$RANGE), vol$vol_func, out$detect.function)$value
}
##########################################################################
###############  Yelloweye rockfish example - Stonewall bank, Oregon coast 2019 ############
if (runyelloweye==TRUE){
  # load yelloweye dataset
  load('data/yelloweye.rda')

  # get the volume estimated
  # for this example, we are removing the seafloor from teh picture.  We know the camera is 0.35 m from teh seafloor, which is tilted up 1.74 degrees
  vol <- get_vol_func(yelloweye$cal,max_extent=7, grid_size=0.05, plotting=TRUE,units='m3', seafloor_position=c(1.74,0,0.35))

  # get density data
  detection_data=prep_detection_data(target_ranges=yelloweye$targets$RANGE,
                                     vol_func=vol$vol_func, nbins=25, method='median', nvals=5, loc_dens=NULL, plotting=TRUE)
  # round inputs to remove glm warning
  detection_data$exp_count <- round(detection_data$exp_count)
  detection_data$obs_count <- round(detection_data$obs_count)

  # fit detection function
  out <- fit_density_function(detection_data,method='logistic gam',formula=NULL, plotting=TRUE)

  # integrate
  eff_vol_yelloweye=integrate(eff_vol_func,lower =0,
                            upper =max(yelloweye$targets$RANGE), vol$vol_func, out$detect.function)$value

}
##########################################################################
###############  Walleye polock example - Gulf of Alaska 2023 ############
if (runpollock==TRUE){
# load pollock dataset
load('data/pollock.rda')

# get the volume estimated
vol <- get_vol_func(pollock$cal,max_extent=8, grid_size=0.05, plotting=TRUE)

# get density data
detection_data=prep_detection_data(target_ranges=pollock$targets$RANGE,
                                   vol_func=vol$vol_func, nbins=25, method='median', nvals=5, loc_dens=NULL, plotting=TRUE)
# round inputs to remove glm warning
detection_data$exp_count <- round(detection_data$exp_count)
detection_data$obs_count <- round(detection_data$obs_count)

# fit detection function
out <- fit_density_function(detection_data,method='logistic gam',formula=NULL, plotting=TRUE)

# integrate
eff_vol_pollock=integrate(eff_vol_func,lower =0,
                            upper =max(pollock$targets$RANGE), vol$vol_func, out$detect.function)$value

}
##################################################################################
####################### Krill example GOA 2015 ###############################
if (runkrill==TRUE){
  # load pollock dataset
  load('data/krill.rda')

  # get the volume estimated
  vol <- get_vol_func(krill$cal,max_extent=20, grid_size=0.1, plotting=TRUE, units='l')

  # get density data
  detection_data=prep_detection_data(target_ranges=krill$targets$Range,
                                     vol_func=vol$vol_func, nbins=25, method='median', nvals=5, loc_dens=NULL, plotting=TRUE)
  # round inputs to remove glm warning
  detection_data$exp_count <- round(detection_data$exp_count)
  detection_data$obs_count <- round(detection_data$obs_count)

  # fit detection function
  out <- fit_density_function(detection_data,method='logistic gam',formula=NULL, plotting=TRUE)

  # integrate
  eff_vol_krill=integrate(eff_vol_func,lower =0,
                            upper =max(krill$targets$Range), vol$vol_func, out$detect.function)$value

}
