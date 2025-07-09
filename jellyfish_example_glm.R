# this is an example script with all the components in a single run, with the goal of
# providing a baseline for comparing a more structured approach - without plotting!

rm(list = ls())
library(StereoCamVolume)

# this data set comes with the package now - replaced with a read function bc I don't know how this runs :-(
load('data/jellyfish.rda')

##########################################################################
# get the volume estimated
vol_func <- get_vol_func(jellyfish$cal,max_extent=7, grid_size=0.05)
# plotting
plot(vol_func$range_centers,vol_func$vol)
x=seq(0, max(vol_func$range_centers),length.out=100)
y=vol_func$vol_func(x)
lines(x,y,col="red")

# try fitting all at once
out1 <- fit_single_logistic(target_ranges=jellyfish$targets$RANGE,vol_func=vol_func$vol_func,nbins=25,plotting=FALSE)
exp_dens=out1$model_output$par[1]
out2 <- fit_single_logistic_glm(target_ranges=jellyfish$targets$RANGE,vol_func=vol_func$vol_func,nbins=25,plotting=FALSE, exp_dens=exp_dens)

exp_dens=out1$model_output$par[1]
exp_count=exp_dens*out1$density_data$vol
plot(out1$density_data$rng,out1$density_data$bin_count/exp_count,xlab='range from camera (m)',ylab='observed/expected')

x=seq(0,max(jellyfish$targets$RANGE),length.out=100)
out1$model_output$par[1]=1
y1=func_single_logistic(out1$model_output$par,x)
L50=-out2$model_output$coefficients[1]/out2$model_output$coefficients[2]
SR=2*log(3)/out2$model_output$coefficients[2]
par_glm=c(1,L50,SR)
y1=func_single_logistic(out1$model_output$par,x)
y2=func_single_logistic(par_glm,x)
lines(x,y1,col="red")
lines(x,y2,col="blue")




