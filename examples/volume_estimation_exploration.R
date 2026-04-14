# this example follows the volume estimation process step by step,including 3D plots for ilustrative purposes

rm(list = ls())

library(StereoCamVolume)
# read in a calibration
Cal=read_calibration('accessory_functions/example_calibrations/matlab_caltech_example2.mat','matlab_caltech')

volume_out <- get_vol_func(jellyfish$cal,max_extent=4, grid_size=0.1, force_origin=FALSE,plotting=FALSE, units='m3')

x=seq(0,0.4,length.out=100)
y1=volume_out$vol_func(x)

plot(x, y1, type = "l", lty = 3)

volume_out <- get_vol_func(jellyfish$cal,max_extent=4, grid_size=0.1, force_origin=TRUE,plotting=FALSE, units='m3')


y2=volume_out$vol_func(x)


abline(h=0)
abline(v=volume_out$x0)
lines(x,y2,col="red")
