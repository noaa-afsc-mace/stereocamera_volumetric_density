# this is an example script with all the components in a single run, with the goal of
# providing a baseline for comparing a more structured approach - without plotting!

rm(list = ls())
library(StereoCamVolume)
library(ggplot2)

# this data set comes with the package now - replaced with a read function bc I don't know how this runs :-(
load('data/jellyfish.rda')

##########################################################################
# get the volume estimated
volume_out <- get_vol_func(jellyfish$cal,max_extent=4, grid_size=0.05, plotting=FALSE)

# run the boostrap
bootstrap_out=bootstrap_effective_volume(target_ranges=jellyfish$targets$RANGE,
                           vol_func=volume_out$vol_func,
                           nbins=15,
                           method='mean',
                           nvals=3, loc_dens=NULL,
                           model_method='logistic glm',
                           formula=cbind(obs_count, exp_count - obs_count) ~ range,
                           boot_n=500,
                           plotting=FALSE)#

# plot histogram of effective volumes

plotdf=data.frame(eff_vol=bootstrap_out$eff_vol)
p1=ggplot(plotdf, aes(x=eff_vol)) +
  geom_histogram(fill="gray45") +
  labs(title = "", x = expression("effective volume (m"^"3)"), y = "frequency")  +
  theme_bw() +
  geom_vline(xintercept=mean(bootstrap_out$eff_vol))+
  xlim(0.11, .23)+
  scale_x_continuous(breaks=c(0.12,0.14,0.16,0.18,0.2,0.22), labels=c(0.12,0.14,0.16,0.18,0.20,0.22))


print(p1)
# plot the lines
# get basic data for plotting
prep_out=prep_detection_data(target_ranges=jellyfish$targets$RANGE,
                             vol_func=volume_out$vol_func, nbins=15, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)

plotdf=data.frame(cbind(x=prep_out$data$range, y=bootstrap_out$fit))
plotdf2=data.frame(cbind(x=prep_out$data$range, y1=bootstrap_out$lower_CI, y2=bootstrap_out$upper_CI))
p2=ggplot(prep_out$data, aes(x=range, y=obs_count/exp_count)) +
  geom_point() +
  geom_line(data=plotdf, mapping=aes(x=x, y = y), color = "black", linewidth = 0.5) +
  geom_line(data=plotdf2, mapping=aes(x=x, y = y1), color = "black", linewidth = 0.5,linetype = "dashed") +
  geom_line(data=plotdf2, mapping=aes(x=x, y = y2), color = "black", linewidth = 0.5,linetype = "dashed") +
  labs(title = "", x = "range from camera (m)", y = "detection probability")  +
  theme_bw()
print(p2)
