# this is an example script with all the components in a single run, with the goal of
# providing a baseline for comparing a more structured approach - without plotting!



rm(list = ls())
library(likelihoodExplore)
library(tidyverse)
library(optimx)


##########################################################################
# get the volume estimated
library(R.matlab)
library(ptinpoly)
library(StereoCamVolume)

# reading in matlab file with stereo calibration parameters
matcal <- readMat("volume_estimation//minicam_09262023_cal_sebastes.mat")
vol_func <- get_vol_func(matcal)

######################################################################
#  TARGET DETECTION FUNCTION ############

# read in fish data
targets<-read.csv("detection_function/targets.csv",header=TRUE)
targets=targets %>% filter(SPECIES_GROUP=="Adult_pollock")
targets$RANGE=targets$RANGE/100# change units from cm to m


dens <- get_dens(targets=targets, vol_func=vol_func$vol_func)

# minimization
out=opm(par=c(500,4,-2), fn=detect_single_logistic, gr = NULL, vol=dens$vol,
        rng=dens$bin_mids,fishcount=dens$bin_counts,method="Nelder-Mead")


# now we multiply the two functions and integrate
# we compute the effective area based on volume and detection functions
params=c(vol_func$p_vol,out$p1,out$p2)
# the parameters here are essentially L50, and SR equivalents, rather than the
# typical slope intercept of the logistic regression.  The third parameter in the optimization is the scale paramter,
# we ignore it so that the function maximizes at 1
# define the complete function - volume 2nd order polynomial and the logistic detection function
eff_vol_func<- function(x,params){(params[1]+params[2]*x+params[3]*x^2)*(1./(1+9^((params[4]-x)/params[5])))}

# integrate
eff_vol=integrate(eff_vol_func,params,lower =0, upper =15)

# now we can turn the counts-by-frame into density-by-frame
count_targets=targets %>% group_by(FRAME_NUMBER) %>% summarise(COUNT=n()) %>% mutate(DENSITY=COUNT/eff_vol$value)
# need to adjust this for "empty" frames
all_frames=seq(min(count_targets$FRAME_NUMBER),max(count_targets$FRAME_NUMBER))
empty_frames=data.frame("FRAME_NUMBER"=setdiff(all_frames,count_targets$FRAME_NUMBER))
empty_frames$COUNT=0
empty_frames$DENSITY=0
count_targets=rbind(count_targets,empty_frames)
summary(count_targets)
print(sum(count_targets$DENSITY))
# 4.462108 should match if nothing changed
