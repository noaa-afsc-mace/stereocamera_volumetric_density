# this is an example script with all the components in a single run, with the goal of
# providing a baseline for comparing a more structured approach - without plotting!



rm(list = ls())

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
targets=targets |> subset(SPECIES_GROUP=="Adult_pollock")
targets$RANGE=targets$RANGE/100# change units from cm to m


dens <- get_dens(targets=targets, vol_func=vol_func$vol_func)

# minimization
out <- optim(par=c(500,4,-2), fn = detect_single_logistic, vol=dens$vol,
              rng=dens$bin_mids,fishcount=dens$bin_counts)

# now we multiply the two functions and integrate
# we compute the effective area based on volume and detection functions
params=c(vol_func$p_vol,out$par[1],out$par[2])
# the parameters here are essentially L50, and SR equivalents, rather than the
# typical slope intercept of the logistic regression.  The third parameter in the optimization is the scale paramter,
# we ignore it so that the function maximizes at 1
# define the complete function - volume 2nd order polynomial and the logistic detection function
eff_vol_func<- function(x,params){
  vol_func$vol_func(x)*(1./(1+9^((params[4]-x)/params[5])))
}

xseq <- seq(0,15, len=1000)
plot(xseq, (params[1]+params[2]*xseq+params[3]*xseq^2))
plot(xseq, (1./(1+9^((params[4]-xseq)/params[5]))))
plot(xseq, eff_vol_func(xseq, params))

# integrate
eff_vol=integrate(eff_vol_func,params,lower =0, upper =15)

# now we can turn the counts-by-frame into density-by-frame
count_targets <- aggregate(FRAME_NUMBER ~ FRAME_NUMBER, data = targets, FUN = length)
count_targets <- as.data.frame(table(targets$FRAME_NUMBER), stringsAsFactors = FALSE)
names(count_targets) <- c("FRAME_NUMBER", "COUNT")
count_targets$FRAME_NUMBER <- as.numeric(count_targets$FRAME_NUMBER)
count_targets$DENSITY <- with(count_targets, COUNT/eff_vol$value)

# need to adjust this for "empty" frames
all_frames <- with(count_targets, seq(min(FRAME_NUMBER),max(FRAME_NUMBER)))
empty_frames <- data.frame("FRAME_NUMBER"=setdiff(all_frames,count_targets$FRAME_NUMBER), COUNT=0, DENSITY=0)
count_targets <- rbind(count_targets,empty_frames)
summary(count_targets)
print(sum(count_targets$DENSITY))
# 4.462108 should match if nothing changed
