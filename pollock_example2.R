# this is an example script with all the components in a single run, with the goal of
# providing a baseline for comparing a more structured approach - without plotting!

rm(list = ls())
library(StereoCamVolume)

# this data set comes with the package now - replaced with a read function bc I don't know how this runs :-(
load('data/pollock.rda')

##########################################################################
# get the volume estimated
vol_func <- get_vol_func(pollock$cal,max_extent=7, grid_size=0.05)
# plotting
plot(vol_func$range_centers,vol_func$vol)
x=seq(0, max(vol_func$range_centers),length.out=100)
y=vol_func$vol_func(x)
lines(x,y,col="red")

dens <- get_dens(targets=pollock$targets, vol_func=vol_func$vol_func)

plot(dens$bin_mids,dens$dens)
# minimization
# out <- optim(par=c(250,3,-2), fn = eval_single_logistic, vol=dens$vol,
#               rng=dens$bin_mids,fishcount=dens$bin_counts)
# y=func_single_logistic(out$par,x)

# out <- optim(par=c(250,1.5,1,3,-2), fn = eval_double_logistic, vol=dens$vol,
#              rng=dens$bin_mids,fishcount=dens$bin_counts)
# y=func_double_logistic(out$par,x)

out <- optim(par=c(250,2,1), fn = eval_normal, vol=dens$vol,
              rng=dens$bin_mids,fishcount=dens$bin_counts)

y=func_normal(out$par,x)

lines(x,y,col="red")


# now we multiply the two functions and integrate
# we compute the effective area based on volume and detection functions
params=c(vol_func$p_vol,out$par)
params[4]=1
# the parameters here are essentially L50, and SR equivalents, rather than the
# typical slope intercept of the logistic regression.  The third parameter in the optimization is the scale paramter,
# we ignore it so that the function maximizes at 1
# define the complete function - volume 2nd order polynomial and the logistic detection function
eff_vol_func<- function(x,params){
  vol_func$vol_func(x)*func_normal(params[4:6],x)
}

# integrate
eff_vol=integrate(eff_vol_func,params,lower =0, upper =15)$value

# now we can turn the counts-by-frame into density-by-frame
count_targets <- as.data.frame(table(pollock$targets$FRAME_NUMBER),
                               stringsAsFactors = FALSE) |>
  setNames(c("FRAME_NUMBER", "COUNT"))
count_targets$FRAME_NUMBER <- as.numeric(count_targets$FRAME_NUMBER)
count_targets$DENSITY <- with(count_targets, COUNT/eff_vol)

# need to adjust this for "empty" frames
all_frames <- with(count_targets, seq(min(FRAME_NUMBER),max(FRAME_NUMBER)))
empty_frames <- data.frame("FRAME_NUMBER"=setdiff(all_frames,count_targets$FRAME_NUMBER), COUNT=0, DENSITY=0)
count_targets <- rbind(count_targets,empty_frames)
count_targets <- count_targets[order(count_targets$FRAME_NUMBER),]
# checks
summary(count_targets)
print(mean(count_targets$DENSITY))

