# this is an example script with all the components in a single run, with the goal of
# providing a baseline for comparing a more structured approach - without plotting!

rm(list = ls())
library(StereoCamVolume)

# this data set comes with the package now - replaced with a read function bc I don't know how this runs :-(
load('data/jellyfish.rda')

##########################################################################
# get the volume estimated
vol <- get_vol_func(jellyfish$cal,max_extent=7, grid_size=0.05, plotting=TRUE)

# get density data
detection_data=prep_detection_data(target_ranges=jellyfish$targets$RANGE,
    vol_func=vol$vol_func, nbins=20, method='median', nvals=5, loc_dens=NULL, plotting=TRUE)

# fit detection function
out <- fit_density_function(detection_data,method='logistic glm',plotting=TRUE)





