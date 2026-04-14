# this is an example script with all the components in a single run, with the goal of
# providing a baseline for comparing a more structured approach - without plotting!

rm(list = ls())
library(StereoCamVolume)
library(ggplot2)
runjellyfish=TRUE
runyelloweye=FALSE
runpollock=FALSE
runkrill=FALSE


if (runjellyfish==TRUE){
##########################################################################
################## Jellyfish (Sarsia sp.) example - Bering Sea 2024 ######

# this data set comes with the package now - replaced with a read function bc I don't know how this runs :-(
load('data/jellyfish.rda')

# get the volume estimated
volume_out <- get_vol_func(jellyfish$cal,max_extent=4, grid_size=0.05, plotting=FALSE, units='m3')
# volume_plot
x=seq(0,2.5,length.out=100)
y=volume_out$vol_func(x)
plotdf=data.frame(cbind(x,y))
p0=ggplot(plotdf, aes(x=x, y=y)) +

  geom_line(color = "black", linewidth = 0.5) +
  labs(title = "", x = "range from camera (m)", y = "change in volume")  +
  theme_bw()
print(p0)

# get density data
prep_out=prep_detection_data(target_ranges=jellyfish$targets$RANGE,
    vol_func=volume_out$vol_func, nbins=15, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)
# plot 1 - density prep for the detection functin input

p1=ggplot(prep_out$data, aes(x=range, y=dens)) +
  geom_bar(stat = "identity") +
  labs(title = "Sarsia sp.", x = "range from camera (m)", y = expression("est. agg. density (#/m"^"3)"))  +
  theme_bw() + geom_hline(yintercept=prep_out$loc_dens)+
  annotate("text", x = 1.15, y = prep_out$loc_dens+80, label = "Est. Max Dens.")
print(p1)

# fit detection function
detection_out <- fit_density_function(prep_out$data,method='logistic glm',formula=NULL, plotting=FALSE)

x=seq(0,3,length.out=100)
y=detection_out$detect.function(x)
plotdf=data.frame(cbind(x,y))
p2=ggplot(prep_out$data, aes(x=range, y=obs_count/exp_count)) + geom_point() +
  geom_line(data=plotdf, mapping=aes(x=x, y = y), color = "black", linewidth = 0.5) +
 labs(title = "", x = "range from camera (m)", y = "detection probability")  +
  theme_bw()
print(p2)

# compute effective volume
eff_vol_jellyfish=eff_vol_compute(integration_range=c(0,max(jellyfish$targets$RANGE)),volume_out$vol_func, detection_out$detect.function)

# plot the area under curve
x=seq(0,2.5,length.out=100)
y=eff_vol_func(x,volume_out$vol_func, detection_out$detect.function)
plotdf=data.frame(cbind(x,y))
p3=ggplot(plotdf, aes(x=x, y=y)) +
  geom_area( fill="gray45", alpha=0.4) +
  geom_line(color = "black", linewidth = 0.5) +
  labs(title = "", x = "range from camera (m)", y = "Volume x detection")  +
  theme_bw() +
  annotate("text", x = 0.85, y = 0.02, label = expression("Effective Volume = 0.14 m"^"3"))
print(p3)

}


##########################################################################
###############  Yelloweye rockfish example - Stonewall bank, Oregon coast 2019 ############
if (runyelloweye==TRUE){
  # load yelloweye dataset
  load('data/yelloweye.rda')

  # get the volume estimated
  # for this example, we are removing the seafloor from teh picture.  We know the camera is 0.35 m from teh seafloor, which is tilted up 1.74 degrees
  volume_out <- get_vol_func(yelloweye$cal,max_extent=7, grid_size=0.05, plotting=FALSE,units='m3', seafloor_position=c(1.74,0,0.35))

  # get density data
  prep_out=prep_detection_data(target_ranges=yelloweye$targets$RANGE,
                                     vol_func=volume_out$vol_func, nbins=25, method='median', nvals=5, loc_dens=NULL, plotting=FALSE)
  # plot 1 - density prep for the detection functin input

  p1=ggplot(prep_out$data, aes(x=range, y=dens)) +
    geom_bar(stat = "identity") +
    labs(title = "Yelloweye Rockfish", x = "range from camera (m)", y = expression("est. agg. density (#/m"^"3)"))  +
    theme_bw() + geom_hline(yintercept=prep_out$loc_dens)+
    annotate("text", x = 5, y = prep_out$loc_dens+5, label = "Est. Max Dens.")
  print(p1)

  # fit detection function
  detection_out <- fit_density_function(prep_out$data,method='logistic gam',formula=NULL, plotting=FALSE)
  x=seq(min(prep_out$data$range),max(prep_out$data$range),length.out=100)
  y=detection_out$detect.function(x)
  plotdf=data.frame(cbind(x,y))
  p2=ggplot(prep_out$data, aes(x=range, y=obs_count/exp_count)) +
    geom_point() +
    geom_line(data=plotdf, mapping=aes(x=x, y = y), color = "black", linewidth = 0.5) +
    labs(title = "", x = "range from camera (m)", y = "detection probability")  +
    theme_bw()
  print(p2)

  # integrate
  eff_vol_yelloweye=integrate(eff_vol_func,lower =0,
                            upper =max(yelloweye$targets$RANGE), volume_out$vol_func, detection_out$detect.function)$value
  eff_vol_label=as.character(round(eff_vol_yelloweye,2) )
  # plot the area under curve
  x=seq(0,max(prep_out$data$range)+0.5,length.out=100)
  y=eff_vol_func(x,volume_out$vol_func, detection_out$detect.function)
  plotdf=data.frame(cbind(x,y))
  p3=ggplot(plotdf, aes(x=x, y=y)) +
    geom_area( fill="gray45", alpha=0.4) +
    geom_line(color = "black", linewidth = 0.5) +
    labs(title = "", x = "range from camera (m)", y = "Volume x detection")  +
    theme_bw() +
    annotate("text", x = 2.7, y = 0.06, label = expression("Effective Volume = 2.69 m"^"3"))
  print(p3)
}
##########################################################################
###############  Walleye polock example - Gulf of Alaska 2023 ############
if (runpollock==TRUE){
# load pollock dataset
load('data/pollock.rda')

# get the volume estimated
volume_out <- get_vol_func(pollock$cal,max_extent=8, grid_size=0.05, plotting=FALSE)

# get density data
prep_out=prep_detection_data(target_ranges=pollock$targets$RANGE,
                                   vol_func=volume_out$vol_func, nbins=25, method='median', nvals=5, loc_dens=NULL, plotting=FALSE)
p1=ggplot(prep_out$data, aes(x=range, y=dens)) +
  geom_bar(stat = "identity") +
  labs(title = "Walleye Pollock", x = "range from camera (m)", y = expression("est. agg. density (#/m"^"3)"))  +
  theme_bw() + geom_hline(yintercept=prep_out$loc_dens)+
  annotate("text", x = 4, y = prep_out$loc_dens-8, label = "Est. Max Dens.")
print(p1)
# fit detection function
detection_out <- fit_density_function(prep_out$data,method='logistic gam',formula=NULL, plotting=FALSE)

# fit detection function
detection_out <- fit_density_function(prep_out$data,method='logistic gam',formula=NULL, plotting=FALSE)
x=seq(min(prep_out$data$range),max(prep_out$data$range),length.out=100)
y=detection_out$detect.function(x)
plotdf=data.frame(cbind(x,y))
p2=ggplot(prep_out$data, aes(x=range, y=obs_count/exp_count)) +
  geom_point() +
  geom_line(data=plotdf, mapping=aes(x=x, y = y), color = "black", linewidth = 0.5) +
  labs(title = "", x = "range from camera (m)", y = "detection probability")  +
  theme_bw()
print(p2)
# integrate
eff_vol_pollock=integrate(eff_vol_func,lower =0,
                            upper =max(pollock$targets$RANGE), volume_out$vol_func, detection_out$detect.function)$value

x=seq(0,max(prep_out$data$range)+0.5,length.out=100)
y=eff_vol_func(x,volume_out$vol_func, detection_out$detect.function)
plotdf=data.frame(cbind(x,y))
p3=ggplot(plotdf, aes(x=x, y=y)) +
  geom_area( fill="gray45", alpha=0.4) +
  geom_line(color = "black", linewidth = 0.5) +
  labs(title = "", x = "range from camera (m)", y = "Volume x detection")  +
  theme_bw() +
  annotate("text", x = 2.7, y = 0.2, label = expression("Effective Volume = 5.86 m"^"3"))
print(p3)
}
##################################################################################
####################### Krill example GOA 2015 ###############################
if (runkrill==TRUE){
  # load pollock dataset
  load('data/krill.rda')

  # get the volume estimated
  volume_out <- get_vol_func(krill$cal,max_extent=20, grid_size=0.5, plotting=FALSE, units='l')
  # volume_plot
  plotdf0=data.frame(cbind(volume_out$range_centers, volume_out$vol))
  colnames(plotdf0) <- c("range","vol")
  x=seq(0,5,length.out=100)
  y=volume_out$vol_func(x)
  plotdf=data.frame(cbind(x,y))
  p0=ggplot(plotdf0,aes(x=range, y=vol)) +
    geom_point() +
    geom_line(data=plotdf, mapping=aes(x=x, y=y),color = "black", linewidth = 0.5) +
    labs(title = "", x = "range from camera (m)", y = "change in volume")  +
    theme_bw()
  print(p0)
  # get density data
  prep_out=prep_detection_data(target_ranges=krill$targets$Range,
                                     vol_func=volume_out$vol_func, nbins=25, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)
  # plot 1 - density prep for the detection functin input

  p1=ggplot(prep_out$data, aes(x=range, y=dens)) +
    geom_bar(stat = "identity") +
    labs(title = "Euphausia sp.", x = "range from camera (m)", y = expression("est. agg. density (#/m"^"3)"))  +
    theme_bw() + geom_hline(yintercept=prep_out$loc_dens)+
    annotate("text", x = 10, y = prep_out$loc_dens-3, label = "Est. Max Dens.") +
    scale_x_continuous(breaks=c(5,10,15), labels=c(0.5,1,1.5))
  print(p1)
  # fit detection function
  detection_out <- fit_density_function(prep_out$data,method='logistic gam',formula=NULL, plotting=FALSE)

  x=seq(min(prep_out$data$range),max(prep_out$data$range),length.out=100)
  y=detection_out$detect.function(x)
  plotdf=data.frame(cbind(x,y))
  p2=ggplot(prep_out$data, aes(x=range, y=obs_count/exp_count)) +
    geom_point() +
    geom_line(data=plotdf, mapping=aes(x=x, y = y), color = "black", linewidth = 0.5) +
    labs(title = "", x = "range from camera (m)", y = "detection probability")  +
    scale_x_continuous(breaks=c(5,10,15), labels=c(0.5,1,1.5)) +
    theme_bw()
  print(p2)
  # integrate
  eff_vol_krill=integrate(eff_vol_func,lower =0,
                            upper =max(krill$targets$Range), volume_out$vol_func, detection_out$detect.function)$value
  # plot the area under curve
  x=seq(0,max(prep_out$data$range)+0.5,length.out=100)
  y=eff_vol_func(x,volume_out$vol_func, detection_out$detect.function)
  y[which(y<0)]=0
  plotdf=data.frame(cbind(x,y))
  p3=ggplot(plotdf, aes(x=x, y=y)) +
    geom_area( fill="gray45", alpha=0.4) +
    geom_line(color = "black", linewidth = 0.5) +
    labs(title = "", x = "range from camera (m)", y = "Volume x detection")  +
    scale_x_continuous(breaks=c(5,10,15), labels=c(0.5,1,1.5)) +
    theme_bw() +
    annotate("text", x = 6.8, y = 2, label = expression("Effective Volume = 0.12 m"^"3"))
  print(p3)
}
