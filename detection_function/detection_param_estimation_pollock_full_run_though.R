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

# reading in matlab file with stereo calibration parameters
matcal <- readMat("volume_estimation//minicam_09262023_cal_sebastes.mat")


# code below is adopted from the camera calibration toolbox 
#http://robots.stanford.edu/cs223b04/JeanYvesCalib/
# cited as Bouguet, J.Y., 2008. Camera calibration toolbox for Matlab [online]. [Available from http://vision.caltech.edu/bouguetj/calib_doc/index.html (accessed September 2008)].

# this code sets up the view "cones"  
# we use this to figure out which points in space are visible to the camera 
normT=8000 # extent of visual field in mm
im_dim=c(2048,1536)

BASE_left = normT*(t(matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1),3,3)) %*% 
                     t(matrix(c(1/matcal$fc.left[1], 0, 0, 0, 1/matcal$fc.left[1], 0, 0, 0, 1),3,3)) %*% 
                     t(matrix(c(1, 0, -matcal$cc.left[1], 0, 1, -matcal$cc.left[2], 0, 0, 1),3,3)) %*% 
                     t(matrix(c(0, im_dim[1]-1, im_dim[1]-1, 0, 0, 0, 0, im_dim[2]-1, im_dim[2]-1, 0, 1, 1, 1, 1, 1),5,3)))
IP_left  = matrix(rbind(BASE_left,matrix(0,3,5),BASE_left),3,15)

BASE_right = normT*(t(matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1),3,3)) %*% 
                      t(matrix(c(1/matcal$fc.right[1], 0, 0, 0, 1/matcal$fc.right[1], 0, 0, 0, 1),3,3)) %*% 
                      t(matrix(c(1, 0, -matcal$cc.right[1], 0, 1, -matcal$cc.right[2], 0, 0, 1),3,3)) %*% 
                      t(matrix(c(0, im_dim[1]-1, im_dim[1]-1, 0, 0, 0, 0, im_dim[2]-1, im_dim[2]-1, 0, 1, 1, 1, 1, 1),5,3)))
IP_right  = matrix(rbind(BASE_right,matrix(0,3,5),BASE_right),3,15)
IP_right = matcal$R %*% (IP_right - matrix(matcal$T,3, 15))

# some of this is redundant, which causes issues later on, so I'm paring it down here
IP_left=IP_left[,c(1,2,4,4,5,7,7,8,10,10,11,13,1,4,7,10)]
IP_right=IP_right[,c(1,2,4,4,5,7,7,8,10,10,11,13,1,4,7,10)]

#originally this is for a down facing camera, so we flip the z and y axes, and reverse the direction of the z axis.  Also, change from mm to m
xl=IP_left[1,]/1000
yl=IP_left[3,]/1000
zl=-IP_left[2,]/1000

xr=IP_right[1,]/1000
yr=IP_right[3,]/1000
zr=-IP_right[2,]/1000


# generate grid points
sc=0.1
# making sure we envelop the entire view cones with points
xg_vec=seq(floor(min(c(xl,xr))*10)/10,ceiling(max(c(xl,xr))*10)/10,by=sc)
yg_vec=seq(floor(min(c(yl,yr))*10)/10,ceiling(max(c(yl,yr))*10)/10,by=sc)
zg_vec=seq(floor(min(c(zl,zr))*10)/10,ceiling(max(c(zl,zr))*10)/10,by=sc)
grid_frame=expand.grid(xg_vec,yg_vec,zg_vec)
colnames(grid_frame) <- c("x","y","z")
# find ones that are within both cones
pos_l=matrix(c(xl,yl,zl),length(xl),3)
pos_r=cbind(xr,yr,zr)
# left cone
grid_full=matrix(c(grid_frame$x,grid_frame$y,grid_frame$z),length(grid_frame$z),3)
faces=t(matrix(c(1,2,3,1,4,3,4,5,3,5,2,3,1,2,4,4,5,2),3,6))
verts_left=pos_l[c(1,9,2,3,6),]
pts_in_left=pip3d(verts_left,faces,grid_full)
in_left=which(pts_in_left==1)
grid_left=grid_full[in_left,]
# right cone
verts_right=pos_r[c(1,9,2,3,6),]
pts_in_right=pip3d(verts_right,faces,grid_left)
in_right=which(pts_in_right==1)
grid_both=grid_left[in_right,]
# now to figure out the "change in volume" function
range_bins=seq(0.5,7.5,by=1)# these are the boundaries - only going out to 13 m as there are edge effects
grid_pt_ranges=sqrt(grid_both[,1]^2+grid_both[,2]^2+grid_both[,3]^2)# euclidean distance to origin (left camera)
vol=vector(mode="numeric",length=7)
v_pt=sc^3
# The key here is that the bin width is equivalent to the unit, (e.g. 1 m).  
# This gives us the change function.
for (i in 1:7){
  # find points in range interval
  ind=which(grid_pt_ranges>range_bins[i] & grid_pt_ranges<=range_bins[i+1])
  # number of points times volume per point
  vol[i]=length(ind)*v_pt
}
cam_range=1:7

# fit 2 order polynomial to this
p_vol= coef(lm(vol ~cam_range +I(cam_range^2)))
disp_range=seq(0,7,length.out = 100)
y=p_vol[1]+p_vol[2]*disp_range+p_vol[3]*disp_range^2

# create the function object for itegration
vol_func<- function(x,p_vol){p_vol[1]+p_vol[2]*x+p_vol[3]*x^2}


######################################################################
#  TARGET DETECTION FUNCTION ############

# read in fish data
targets<-read.csv("detection_function/targets.csv",header=TRUE)
targets=targets %>% filter(SPECIES_GROUP=="Adult_pollock")
targets$RANGE=targets$RANGE/100# change units from cm to m

# here we bin the data for fitting the model.  This involves estimation the volume 
# for each range bin as well as how many animals are in the bin.
nbins=25
bin_edges=seq(min(targets$RANGE),max(targets$RANGE),length.out=nbins+1)
bin_mids=bin_edges[1:nbins]+(bin_edges[2]-bin_edges[1])/2

vol=numeric(length=length(bin_mids))
bin_counts=numeric(length=length(bin_mids))
for (i in 1:length(bin_mids)){
  # this is integrating how much volume is in each bin to get to density
  v=integrate(vol_func,p_vol,lower = bin_edges[i], upper = bin_edges[i+1])
  vol[i]=v$value
  ind=which(targets$RANGE>bin_edges[i] & targets$RANGE<=bin_edges[i+1])
  bin_counts[i]=length(ind)
}

dens=bin_counts/vol

# logistic
detect_single_logistic = function(parm,vol,rng,fishcount){
  like=numeric(length=length(rng))
  for (i in 1:length(rng)){
    # expected number of fish per m3
    d_hat=parm[1]/(1+9^((parm[2]-rng[i])/parm[3]))
    # poisson likelihood?
    # ll=likpois(round(fishcount[i]/vol[i]), round(d_hat), log = TRUE)*fishcount[i]*-1
    # normal?
    ll=liknorm(fishcount[i]/vol[i], d_hat, log = TRUE)*fishcount[i]*-1 
    # key part here is to multiply the likelihood by the number of observations in the bin.
    # this way bins with more fish in them drive the fit
    like[i]=ll
  }
  return(sum(like))
}
# minimization
out=opm(par=c(500,4,-2), fn=detect_single_logistic, gr = NULL, vol=vol,
        rng=bin_mids,fishcount=bin_counts,method="Nelder-Mead")


# now we multiply the two functions and integrate
# we compute the effective area based on volume and detection functions
params=c(p_vol,out$p1,out$p2)
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
