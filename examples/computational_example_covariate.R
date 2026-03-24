# this is an example script with all the components in a single run, with the goal of
# providing a baseline for comparing a more structured approach - without plotting!

rm(list = ls())
library(StereoCamVolume)
library(ggplot2)

# greenstrip target data
targets<-read.csv("data-raw/greenstriped_targets2.csv",header=TRUE)
targets$range=targets$range/100# change units from cm to m
# greenstripd deployment data
deployments<-read.csv("data-raw/drops_for_2019_greenstripe.csv",header=TRUE)
ind=which(deployments$keep>0)
deployments=deployments[ind,]
# plot raw relationship between sighting distance and averge fish distance
lm_r = lm(formula = mean_fish_dist ~ farthest_dist, data = deployments)
x=c(min(deployments$farthest_dist),max(deployments$farthest_dist))
plotdf=data.frame(x=x,y=predict(lm_r, newdata = data.frame(farthest_dist=x)))
p1=ggplot(deployments, aes(x=farthest_dist, y=mean_fish_dist)) +
  geom_point(shape = 1, size = 3) +
  geom_line(data=plotdf, mapping=aes(x=x, y=y), color = "black", linewidth = 0.5) +
  labs(title = "a", x = "max sighting distance (m)", y = expression("mean target distance (m)"),color="Legend")  +
  theme_bw()
# print(p1)
#density_data=data.frame(range=bin_mids,dens,exp_count,obs_count=bin_counts)
prep_data_frame=data.frame(range = numeric(),
                             dens = numeric(),
                             exp_count = integer(),
                             obs_count = integer(),
                           cov_value = numeric())
for (i in 1:length(deployments$deployment_id)){
  # this next piece has to go in a loop, because we have different trigcams
  # get the volume estimated
  # pick the appropriate cal
  if (deployments$unit[i]==101){
    Cal=readCalibration('C:/Users/kresimir.williams/Work/trigcam/2019_coop_analysis/KW_analysis/unit_1_regcal.mat','matlab_caltech')
  }
  else if (deployments$unit[i]==102){
    Cal=readCalibration('C:/Users/kresimir.williams/Work/trigcam/2019_coop_analysis/KW_analysis/unit_2_regcal.mat','matlab_caltech')
  }
  else if (deployments$unit[i]==103){
    Cal=readCalibration('C:/Users/kresimir.williams/Work/trigcam/2019_coop_analysis/KW_analysis/unit_3_regcal.mat','matlab_caltech')
  }
  else if (deployments$unit[i]==105){
    Cal=readCalibration('C:/Users/kresimir.williams/Work/trigcam/2019_coop_analysis/KW_analysis/unit_6_regcal.mat','matlab_caltech')
  }
  else if (deployments$unit[i]==108){
    Cal=readCalibration('C:/Users/kresimir.williams/Work/trigcam/2019_coop_analysis/KW_analysis/unit_8_regcal.mat','matlab_caltech')
  }
  seafloor_position=c(deployments$angle[i],0,deployments$cam_height[i])
  volume_out <- get_vol_func(Cal,max_extent=5, grid_size=0.05, plotting=FALSE, units='m3', seafloor_position=seafloor_position)
  # get density data
  ind=which(targets$deployment_id==deployments$assembled[i])
  if (length(ind)>0){
    targets_for_dep=targets[ind,]
    prep_data=prep_detection_data(target_ranges=targets_for_dep$range,
                                  vol_func=volume_out$vol_func, nbins=15, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)
    prep_data$data$cov_factor=deployments$farthest_dist[i]
    prep_data_frame=rbind(prep_data_frame,prep_data$data)
    }
  }

# fit detection function
detection_out <- fit_density_function_covariate(prep_data_frame,method='logistic glm',formula=NULL, plotting=FALSE)

cov_levels=c(2,3,4)
x=seq(0,5,length.out=100)

y1=detection_out$detect.function(x,cov_levels[1])
y2=detection_out$detect.function(x,cov_levels[2])
y3=detection_out$detect.function(x,cov_levels[3])
plotdf2=data.frame(cbind(x=x, y1=y1, y2=y2, y3=y3))
p2=ggplot(plotdf2, aes(x=x)) +
  geom_line(data=plotdf2, mapping=aes(y = y1,color = "MSD = 2"),  linewidth = 1) +
  geom_line(data=plotdf2, mapping=aes(y = y2,color = "MSD = 3"),  linewidth = 1) +
  geom_line(data=plotdf2, mapping=aes(y = y3,color = "MSD = 4"),  linewidth = 1) +
  scale_color_manual(name = "", values = c("MSD = 2" = "black", "MSD = 3" = "gray45","MSD = 4" = "gray85")) +
  # annotate("text", x = 2, y = 0.8, label = "Vis cov = 4 m")+
  # annotate("text", x = 2, y = 0.65, label =  "Vis cov = 3 m")+
  # annotate("text", x = 2, y = 0.50, label =  "Vis cov = 2 m")+
  labs(title = "b", x = "range from camera (m)", y = "detection probability")  +
  theme_bw()+
  theme(legend.position="bottom")
# print(p2)

# now a figure that has change in effective volume with covariate
cov_factor=seq(1,4.5,by=0.25)
eff_volume=vector(mode="numeric",length=length(cov_factor))
Cal=readCalibration('C:/Users/kresimir.williams/Work/trigcam/2019_coop_analysis/KW_analysis/unit_1_regcal.mat','matlab_caltech')
seafloor_position=c(mean(deployments$angle),0,mean(deployments$cam_height))
volume_out <- get_vol_func(Cal,max_extent=5, grid_size=0.05, plotting=FALSE, units='m3')

for (i in 1:length(cov_factor)){
  int_obj=integrate(eff_vol_func_covariate,lower =0,upper=5, cov_factor[i], volume_out$vol_func, detection_out$detect.function)
eff_volume[i]=int_obj$value
}
plotdf=data.frame(cov_factor=cov_factor,eff_volume=eff_volume)
p3=ggplot(plotdf, aes(x=cov_factor, y=eff_volume)) +
  geom_point(shape = 1, size = 3) +
  geom_line(data=plotdf, mapping=aes(x=cov_factor, y=eff_volume), color = "black", linewidth = 0.5) +
  labs(title = "c", x = "max sighting distance (m)", y = expression("effective volume (m"^"3"~")"),color="Legend")  +
  theme_bw()
# print(p3)


library(patchwork)
p1 | p2 | p3 # All three plots in one row

