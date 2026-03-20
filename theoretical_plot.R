# this is an example script with all the components in a single run, with the goal of
# providing a baseline for comparing a more structured approach - without plotting!

rm(list = ls())
library(StereoCamVolume)
library(ggplot2)

##########################################################################
################## Jellyfish (Sarsia sp.) example - Bering Sea 2024 ######

# this data set comes with the package now - replaced with a read function bc I don't know how this runs :-(
load('data/jellyfish.rda')

# get the volume estimated
volume_out <- get_vol_func(jellyfish$cal,max_extent=5, grid_size=0.05, plotting=FALSE, units='m3')

# get density data
prep_data=prep_detection_data(target_ranges=jellyfish$targets$RANGE+0.5,
    vol_func=volume_out$vol_func, nbins=8, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)

# fit detection function
detection_out <- fit_density_function(prep_data$data,method='logistic glm',formula=NULL, plotting=FALSE)

# integrate
eff_vol_jellyfish=integrate(eff_vol_func,lower =0,
                            upper =max(jellyfish$targets$RANGE), volume_out$vol_func, detection_out$detect.function)$value
# plot the effective volume function for kicks
cam_range=seq(0,3,length.out=100)
vol=volume_out$vol_func(cam_range)
detection=detection_out$detect.function(cam_range)
eff_vol=eff_vol_func(cam_range,volume_out$vol_func, detection_out$detect.function)
plotdf=data.frame(cbind(cam_range,vol,detection,eff_vol))


p3 = ggplot(plotdf, aes(x = cam_range)) +
  # 1. Map aesthetics to labels
  geom_area(aes(y = eff_vol, fill = "Effective Volume"), alpha = 0.4) +
  geom_line(aes(y = detection, color = "Detection"), linewidth = 0.5) +
  geom_line(aes(y = vol/5, color = "Volume Change"), linewidth = 1, linetype = "dashed") +

  # 2. Set colors and use the SAME name for both scales
  scale_fill_manual(name = "", values = c("Effective Volume" = "gray45")) +
  scale_color_manual(name = "", values = c("Detection" = "black", "Volume Change" = "gray30")) +

  # 3. Formatting
  labs(x = "range from camera (m)", y = "Detection probability") +
  scale_y_continuous(
    breaks = c(0, 0.5, 1.0),
    labels = c(0, 0.5, 1.0),
    sec.axis = sec_axis(~.*5, name = expression("Change in Volume (m"^3*")"))
  ) +
  theme_bw() +
  theme(legend.position="bottom")

print(p3)

