# this is an example script with all the components in a single run, with the goal of
# providing a baseline for comparing a more structured approach - without plotting!

rm(list = ls())
library(StereoCamVolume)
library(ggplot2)
library(patchwork)
runjellyfish=TRUE
runyelloweye=TRUE
runpollock=TRUE
runkrill=TRUE


if (runjellyfish==TRUE){
##########################################################################
################## Jellyfish (Sarsia sp.) example - Bering Sea 2024 ######

# this data set comes with the package now - replaced with a read function bc I don't know how this runs :-(
load('data/jellyfish.rda')

# get the volume estimated
volume_out <- get_vol_func(jellyfish$cal,max_extent=4, grid_size=0.05, plotting=FALSE, units='m3')
# volume_plot
# x=seq(0,2.5,length.out=100)
# y=volume_out$vol_func(x)
# plotdf=data.frame(cbind(x,y))
# p0=ggplot(plotdf, aes(x=x, y=y)) +
#
#   geom_line(color = "black", linewidth = 0.5) +
#   labs(title = "", x = "range from camera (m)", y = "change in volume")  +
#   theme_bw()
# print(p0)

# get density data
prep_out=prep_detection_data(target_ranges=jellyfish$targets$RANGE,
    vol_func=volume_out$vol_func, nbins=15, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)
# plot 1 - density prep for the detection functin input

p1=ggplot(prep_out$data, aes(x=range, y=dens)) +
  geom_bar(stat = "identity") +
  labs(title = "Sarsia sp.", x = "range from camera (m)", y = expression("est. agg. density (#/m"^"3)"))  +
  theme_bw() + geom_hline(yintercept=prep_out$loc_dens)+
  annotate("text", x = 1.15, y = prep_out$loc_dens+80, label = "Est. Max Dens.")
#print(p1)

# fit detection function
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
eff_vol_mean=mean(bootstrap_out$eff_vol)
eff_vol_CV=sd(bootstrap_out$eff_vol)/eff_vol_mean*100
p3=ggplot(plotdf, aes(x=eff_vol)) +
  geom_histogram(fill="gray45") +
  labs(title = "", x = expression("effective volume (m"^"3"~")"), y = "frequency")  +
  theme_bw() +
  geom_vline(xintercept=mean(bootstrap_out$eff_vol))+
  annotate("text", x = 0.16, y = 5, label = paste("Mean Eff. Vol = ",round(mean(bootstrap_out$eff_vol),3)))+
  xlim(0.11, .23)+
  scale_x_continuous(breaks=c(0.12,0.14,0.16,0.18,0.2,0.22), labels=c(0.12,0.14,0.16,0.18,0.20,0.22))


#print(p1)
# plot the lines
# get basic data for plotting
#prep_out=prep_detection_data(target_ranges=jellyfish$targets$RANGE,
 #                            vol_func=volume_out$vol_func, nbins=15, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)

plotdf=data.frame(cbind(x=prep_out$data$range, y=bootstrap_out$fit))
plotdf2=data.frame(cbind(x=prep_out$data$range, y1=bootstrap_out$lower_CI, y2=bootstrap_out$upper_CI))
p2=ggplot(prep_out$data, aes(x=range, y=obs_count/exp_count)) +
  geom_point() +
  geom_line(data=plotdf, mapping=aes(x=x, y = y, color = "GLM fit"), linewidth = 0.5) +
  geom_line(data=plotdf2, mapping=aes(x=x, y = y1, color="95 % CI"), linewidth = 0.5,linetype = "dotted") +
  geom_line(data=plotdf2, mapping=aes(x=x, y = y2), color = "black", linewidth = 0.5,linetype = "dotted") +
  labs(title = "", x = "range from camera (m)", y = "detection probability")  +
  scale_color_manual(name = "", values = c("GLM fit" = "black","95 % CI"="black"))+
  theme_bw()+

  theme(legend.position = "inside",legend.position.inside=c(1,1),legend.justification = c(1.05, 1.05))

#print(p2)

}


##########################################################################
###############  Yelloweye rockfish example - Stonewall bank, Oregon coast 2019 ############
if (runyelloweye==TRUE){
  # load yelloweye dataset
  load('data/yelloweye.rda')

  # get the volume estimated
  # for this example, we are removing the seafloor from the picture.  We know the camera is 0.35 m from teh seafloor, which is tilted up 1.74 degrees
  volume_out <- get_vol_func(yelloweye$cal,max_extent=7, grid_size=0.05, plotting=FALSE,units='m3', seafloor_position=c(1.74,0,0.35))

  prep_out=prep_detection_data(target_ranges=yelloweye$targets$RANGE,
                               vol_func=volume_out$vol_func, nbins=25, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)
  # plot 1 - density prep for the detection function input

  p4=ggplot(prep_out$data, aes(x=range, y=dens)) +
    geom_bar(stat = "identity") +
    labs(title = "Yelloweye Rockfish", x = "range from camera (m)", y = expression("est. agg. density (#/m"^"3)"))  +
    theme_bw() + geom_hline(yintercept=prep_out$loc_dens)+
    annotate("text", x = 5, y = prep_out$loc_dens+5, label = "Est. Max Dens.")
  # print(p1)
  # run the boostrapn bins=25, method='median', nvals=5,
  bootstrap_out=bootstrap_effective_volume(target_ranges=yelloweye$targets$RANGE,
                                           vol_func=volume_out$vol_func,
                                           nbins=25,
                                           method='mean',
                                           nvals=3, loc_dens=NULL,
                                           model_method='logistic gam',
                                           formula=cbind(obs_count, exp_count - obs_count) ~ range,
                                           boot_n=500,
                                           plotting=FALSE)#

  # plot histogram of effective volumes

  plotdf=data.frame(eff_vol=bootstrap_out$eff_vol)
  eff_vol_mean=mean(bootstrap_out$eff_vol)
  eff_vol_CV=sd(bootstrap_out$eff_vol)/eff_vol_mean*100
  p6=ggplot(plotdf, aes(x=eff_vol)) +
    geom_histogram(fill="gray45") +
    labs(title = "", x = expression("effective volume (m"^"3"~")"), y = "frequency") +
    theme_bw() +
    geom_vline(xintercept=mean(bootstrap_out$eff_vol))+
    annotate("text", x = mean(bootstrap_out$eff_vol), y = 5, label = paste("Mean Eff. Vol = ",round(mean(bootstrap_out$eff_vol),3)))
    # xlim(0.11, .23)
    #scale_x_continuous(breaks=c(0.12,0.14,0.16,0.18,0.2,0.22), labels=c(0.12,0.14,0.16,0.18,0.20,0.22))


  # print(p1)
  # plot the lines
  # get basic data for plotting
  prep_out=prep_detection_data(target_ranges=yelloweye$targets$RANGE,
                               vol_func=volume_out$vol_func, nbins=25, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)

  plotdf=data.frame(cbind(x=prep_out$data$range, y=bootstrap_out$fit))
  plotdf2=data.frame(cbind(x=prep_out$data$range, y1=bootstrap_out$lower_CI, y2=bootstrap_out$upper_CI))
  p5=ggplot(prep_out$data, aes(x=range, y=obs_count/exp_count)) +
    geom_point() +
    geom_line(data=plotdf, mapping=aes(x=x, y = y, color = "GAM fit"), linewidth = 0.5) +
    geom_line(data=plotdf2, mapping=aes(x=x, y = y1, color="95 % CI"), linewidth = 0.5,linetype = "dotted") +
    geom_line(data=plotdf2, mapping=aes(x=x, y = y2), color = "black", linewidth = 0.5,linetype = "dotted") +
    labs(title = "", x = "range from camera (m)", y = "detection probability")  +
    scale_color_manual(name = "", values = c("GAM fit" = "black","95 % CI"="black"))+
    theme_bw()+

    theme(legend.position = "inside",legend.position.inside=c(1,1),legend.justification = c(1.05, 1.05))

   #print(p4)
   # print(p5)
   # print(p6)

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
p7=ggplot(prep_out$data, aes(x=range, y=dens)) +
  geom_bar(stat = "identity") +
  labs(title = "Walleye Pollock", x = "range from camera (m)", y = expression("est. agg. density (#/m"^"3)"))  +
  theme_bw() + geom_hline(yintercept=prep_out$loc_dens)+
  annotate("text", x = 4, y = prep_out$loc_dens-8, label = "Est. Max Dens.")

# run the boostrapn
bootstrap_out=bootstrap_effective_volume(target_ranges=pollock$targets$RANGE,
                                         vol_func=volume_out$vol_func,
                                         nbins=25,
                                         method='mean',
                                         nvals=3, loc_dens=NULL,
                                         model_method='logistic gam',
                                         formula=cbind(obs_count, exp_count - obs_count) ~ range,
                                         boot_n=500,
                                         plotting=FALSE)#

# plot histogram of effective volumes

plotdf=data.frame(eff_vol=bootstrap_out$eff_vol)
eff_vol_mean=mean(bootstrap_out$eff_vol)
eff_vol_CV=sd(bootstrap_out$eff_vol)/eff_vol_mean*100
p9=ggplot(plotdf, aes(x=eff_vol)) +
  geom_histogram(fill="gray45") +
  labs(title = "", x = expression("effective volume (m"^"3"~")"), y = "frequency")  +
  theme_bw() +
  geom_vline(xintercept=mean(bootstrap_out$eff_vol))+
  annotate("text", x = mean(bootstrap_out$eff_vol), y = 5, label = paste("Mean Eff. Vol = ",round(mean(bootstrap_out$eff_vol),3)))
# xlim(0.11, .23)
#scale_x_continuous(breaks=c(0.12,0.14,0.16,0.18,0.2,0.22), labels=c(0.12,0.14,0.16,0.18,0.20,0.22))


#print(p1)
# plot the lines
# get basic data for plotting
prep_out=prep_detection_data(target_ranges=pollock$targets$RANGE,
                             vol_func=volume_out$vol_func, nbins=25, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)

plotdf=data.frame(cbind(x=prep_out$data$range, y=bootstrap_out$fit))
plotdf2=data.frame(cbind(x=prep_out$data$range, y1=bootstrap_out$lower_CI, y2=bootstrap_out$upper_CI))
p8=ggplot(prep_out$data, aes(x=range, y=obs_count/exp_count)) +
  geom_point() +
  geom_line(data=plotdf, mapping=aes(x=x, y = y, color = "GAM fit"), linewidth = 0.5) +
  geom_line(data=plotdf2, mapping=aes(x=x, y = y1, color="95 % CI"), linewidth = 0.5,linetype = "dotted") +
  geom_line(data=plotdf2, mapping=aes(x=x, y = y2), color = "black", linewidth = 0.5,linetype = "dotted") +
  labs(title = "", x = "range from camera (m)", y = "detection probability")  +
  scale_color_manual(name = "", values = c("GAM fit" = "black","95 % CI"="black"))+
  theme_bw()+

  theme(legend.position = "inside",legend.position.inside=c(1,1),legend.justification = c(1.05, 1.05))

#print(p2)

}
##################################################################################
####################### Krill example GOA 2015 ###############################
if (runkrill==TRUE){
  # load pollock dataset
  load('data/krill.rda')

  # get the volume estimated
  volume_out <- get_vol_func(krill$cal,max_extent=20, grid_size=0.5, plotting=FALSE, units='l')
  prep_out=prep_detection_data(target_ranges=krill$targets$Range,
                               vol_func=volume_out$vol_func, nbins=25, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)
  # plot 1 - density prep for the detection functin input

  p10=ggplot(prep_out$data, aes(x=range, y=dens)) +
    geom_bar(stat = "identity") +
    labs(title = "Euphausia sp.", x = "range from camera (m)", y = expression("est. agg. density (#/m"^"3)"))  +
    theme_bw() + geom_hline(yintercept=prep_out$loc_dens)+
    annotate("text", x = 10, y = prep_out$loc_dens-3, label = "Est. Max Dens.") +
    scale_x_continuous(breaks=c(5,10,15), labels=c(0.5,1,1.5))
  # run the boostrapn
  bootstrap_out=bootstrap_effective_volume(target_ranges=krill$targets$Range,
                                           vol_func=volume_out$vol_func,
                                           nbins=25,
                                           method='mean',
                                           nvals=3, loc_dens=NULL,
                                           model_method='logistic glm',
                                           formula=cbind(obs_count, exp_count - obs_count) ~ range,
                                           boot_n=500,
                                           plotting=FALSE)#

  # plot histogram of effective volumes

  plotdf=data.frame(eff_vol=bootstrap_out$eff_vol)
  eff_vol_mean=mean(bootstrap_out$eff_vol)
  eff_vol_CV=sd(bootstrap_out$eff_vol)/eff_vol_mean*100
  p12=ggplot(plotdf, aes(x=eff_vol)) +
    geom_histogram(fill="gray45") +
    labs(title = "", x = expression("effective volume (m"^"3"~")"), y = "frequency")  +
    theme_bw() +
    geom_vline(xintercept=mean(bootstrap_out$eff_vol))+
    annotate("text", x = mean(bootstrap_out$eff_vol), y = 5, label = paste("Mean Eff. Vol = ",round(mean(bootstrap_out$eff_vol/1000),3)))+
  # xlim(0.11, .23)
   scale_x_continuous(breaks=c(110,120,130,140), labels=c(0.11,0.12,0.13,0.14))


  #print(p1)
  # plot the lines
  # get basic data for plotting
  prep_out=prep_detection_data(target_ranges=krill$targets$Range,
                               vol_func=volume_out$vol_func, nbins=25, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)

  plotdf=data.frame(cbind(x=prep_out$data$range, y=bootstrap_out$fit))
  plotdf2=data.frame(cbind(x=prep_out$data$range, y1=bootstrap_out$lower_CI, y2=bootstrap_out$upper_CI))
  p11=ggplot(prep_out$data, aes(x=range, y=obs_count/exp_count)) +
    geom_point() +
    geom_line(data=plotdf, mapping=aes(x=x, y = y, color = "GLM fit"), linewidth = 0.5) +
    geom_line(data=plotdf2, mapping=aes(x=x, y = y1, color="95 % CI"), linewidth = 0.5,linetype = "dotted") +
    geom_line(data=plotdf2, mapping=aes(x=x, y = y2), color = "black", linewidth = 0.5,linetype = "dotted") +
    labs(title = "", x = "range from camera (m)", y = "detection probability")  +
    scale_color_manual(name = "", values = c("GLM fit" = "black","95 % CI"="black"))+
    theme_bw()+

    theme(legend.position = "inside",legend.position.inside=c(1,1),legend.justification = c(1.05, 1.05))

  #print(p2)


}

(p1|p2|p3) / (p4|p5|p6) / (p7|p8|p9) / (p10|p11|p12)
