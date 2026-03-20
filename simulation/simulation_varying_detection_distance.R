rm(list = ls())
library(plot3D)
library(ptinpoly)
library(ggplot2)


load('data/pollock.rda')

# get the volume estimated
volume_out <- get_vol_func(pollock$cal,max_extent=7, grid_size=0.05, plotting=FALSE)

# define detection probability parameters
d50=3# this is the range at which 50% of criters are detection

sl=-1# this is the slope of detection curve, equivalent to selection range from trawl selectivity

scaling=1000 # appropriate for target positions in m
# loop starts here
n = 100000 # number of targets - in this case its targets over a set of frames, so much denser than what you would expect
reps=100
#########################################################################

bins=20
d50=2#
estimated_density1=numeric(length=reps)

for (i in 1:reps){

  # generate density data
  # now we add some random targets to the general space
  xf=runif(n)*10-5
  yf=runif(n)*10
  zf=runif(n)*10-5
  ###### check random 3d distribution ###############

  # find ones that are in both cones
  target_positions=matrix(c(xf,yf,zf),length(xf),3)
  targets_in_both=find_targets_in_view(pollock$cal,scaling, target_positions)

  # estimate probability of detection
  fish_range=sqrt(targets_in_both[,1]^2+targets_in_both[,2]^2+targets_in_both[,3]^2)
  prob=1./(1+9^((d50-fish_range)/sl))
  detected=numeric(length(prob))
  for (j in 1:length(prob)){detected[j]=rbinom(1,1,prob[j])}
  detected_fish_ranges=fish_range[which(detected==1)]
  # run our package estimation
  prep_out=prep_detection_data(target_ranges=detected_fish_ranges,
                                     vol_func=volume_out$vol_func, nbins=bins, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)

  # fit detection function
  detection_out <- fit_density_function(prep_out$data,method='logistic glm',formula=NULL, plotting=FALSE)

  # integrate
  eff_vol_sim=integrate(eff_vol_func,lower =0,
                              upper =max(fish_range), volume_out$vol_func, detection_out$detect.function)$value
  # estimate density
  estimated_density1[i]=length(detected_fish_ranges)/eff_vol_sim
}
#########################################################################

bins=15
d50=3#
estimated_density2=numeric(length=reps)

for (i in 1:reps){

  # generate density data
  # now we add some random targets to the general space
  xf=runif(n)*10-5
  yf=runif(n)*10
  zf=runif(n)*10-5
  ###### check random 3d distribution ###############

  # find ones that are in both cones
  target_positions=matrix(c(xf,yf,zf),length(xf),3)
  targets_in_both=find_targets_in_view(pollock$cal,scaling, target_positions)

  # estimate probability of detection
  fish_range=sqrt(targets_in_both[,1]^2+targets_in_both[,2]^2+targets_in_both[,3]^2)
  prob=1./(1+9^((d50-fish_range)/sl))
  detected=numeric(length(prob))
  for (j in 1:length(prob)){detected[j]=rbinom(1,1,prob[j])}
  detected_fish_ranges=fish_range[which(detected==1)]
  # run our package estimation
  prep_out=prep_detection_data(target_ranges=detected_fish_ranges,
                                     vol_func=volume_out$vol_func, nbins=bins, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)

  # fit detection function
  detection_out <- fit_density_function(prep_out$data,method='logistic glm',formula=NULL, plotting=FALSE)


  # integrate
  eff_vol_sim=integrate(eff_vol_func,lower =0,
                        upper =max(fish_range), volume_out$vol_func, detection_out$detect.function)$value
  # estimate density
  estimated_density2[i]=length(detected_fish_ranges)/eff_vol_sim
}
#########################################################################

bins=15
d50=4
estimated_density3=numeric(length=reps)

for (i in 1:reps){

  # generate density data
  # now we add some random targets to the general space
  xf=runif(n)*10-5
  yf=runif(n)*10
  zf=runif(n)*10-5
  ###### check random 3d distribution ###############

  # find ones that are in both cones
  target_positions=matrix(c(xf,yf,zf),length(xf),3)
  targets_in_both=find_targets_in_view(pollock$cal,scaling, target_positions)

  # estimate probability of detection
  fish_range=sqrt(targets_in_both[,1]^2+targets_in_both[,2]^2+targets_in_both[,3]^2)
  prob=1./(1+9^((d50-fish_range)/sl))
  detected=numeric(length(prob))
  for (j in 1:length(prob)){detected[j]=rbinom(1,1,prob[j])}
  detected_fish_ranges=fish_range[which(detected==1)]
  # run our package estimation
  prep_out=prep_detection_data(target_ranges=detected_fish_ranges,
                                     vol_func=volume_out$vol_func, nbins=bins, method='mean', nvals=3, loc_dens=NULL, plotting=FALSE)

  # fit detection function
  detection_out <- fit_density_function(prep_out$data,method='logistic glm',formula=NULL, plotting=FALSE)


  # integrate
  eff_vol_sim=integrate(eff_vol_func,lower =0,
                        upper =max(fish_range), volume_out$vol_func, detection_out$detect.function)$value
  # estimate density
  estimated_density3[i]=length(detected_fish_ranges)/eff_vol_sim
}

df=data.frame(rbind(cbind(rep(1,length(estimated_density1)),estimated_density1),cbind(rep(2,length(estimated_density2)),estimated_density2),cbind(rep(3,length(estimated_density3)),estimated_density3)))
colnames(df) <- c("bin","val")
df$val=df$val/100
p1 <- ggplot(data=df, aes(group=bin, y=val)) +
  geom_boxplot()+labs(title = "", x = "bin size", y = "est. dens. / true dens")  +
 theme_bw() +scale_x_discrete(expand = c(0.1, 0),name="detection distance",
            limits=c(-0.25,0,.25),labels=c("2 m","3 m","4 m")) + geom_hline(yintercept=1)
print(p1)
