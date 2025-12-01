# this is an example script with all the components in a single run, with the goal of
# providing a baseline for comparing a more structured approach - without plotting!

rm(list = ls())
library(StereoCamVolume)

# this data set comes with the package now - replaced with a read function bc I don't know how this runs :-(
load('data/jellyfish.rda')

##########################################################################
# get the volume estimated
vol <- get_vol_func(jellyfish$cal,max_extent=7, grid_size=0.05, plotting=FALSE)

# run the boostrap
eff_vol=bootstrap_effective_volume(target_ranges=jellyfish$targets$RANGE,
                           vol_func=vol$vol_func,
                           nbins=15,
                           method='median',
                           nvals=5, loc_dens=NULL,
                           model_method='logistic glm',
                           formula=cbind(obs_count, exp_count - obs_count) ~ range,
                           boot_n=500,
                           plotting=TRUE)#cbind(obs_count, exp_count - obs_count) ~ range

# plot histogram of effective volumes
hist(eff_vol,xlab='effective volume (m^3)',ylab='frequency')
abline(v=mean(eff_vol),col="red", lty=2)
