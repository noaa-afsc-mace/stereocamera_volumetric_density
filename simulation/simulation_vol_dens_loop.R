rm(list = ls())
library(plot3D)
library(ptinpoly)

# these are the vertices of our view cones. "l" is for left, "r" for right
xl=c(0,-1,1,0,0,1, 1,0,0,1, -1,0,0,-1,-1,0)*10
zl=c(0,1,1,0,0,1, -1,0,0,-1, -1,0,0,-1,1,0)*10
yl=c(0,1.5,1.5,0,0,1.5,1.5,0,0,1.5,1.5,0,0,1.5,1.5,0)*10

# right view cone is the same with offset in x.  A more realistic scenario would be to use a transformation 
# and rotation matrix from a calibration
xr=xl+1
yr=yl
zr=zl

# make a nice 3d plot
lines3D(xl,yl,zl,xlab = "x (m)", ylab = "y (m)", zlab = "z (m)", col="red",ticktype = "detailed")
lines3D(xr,yr,zr, col="green",add = TRUE)

# for volume comp, make a 3D grid of points
# grid point unit scale, in m
sc=0.5
xg_vec=seq(-11,11,by=sc)
yg_vec=seq(0,15,by=sc)
zg_vec=seq(-11,11,by=sc)
grid_frame=expand.grid(xg_vec,yg_vec,zg_vec)
colnames(grid_frame) <- c("x","y","z")
# plot these if you want, but there are many!
#points3D(grid_frame$x,grid_frame$y,grid_frame$z,col="blue")

# each point is this much volume
v_pt=sc^3
# organize stuff into matrices for "in hull" operator
# for more 
pos_l=matrix(c(xl,yl,zl),length(xl),3)
pos_r=cbind(xr,yr,zr)
grid_full=matrix(c(grid_frame$x,grid_frame$y,grid_frame$z),length(grid_frame$z),3)
faces=t(matrix(c(1,2,3,1,3,4,1,4,5,1,5,2,2,3,5,3,4,5),3,6))
verts_left=pos_l[c(1,2,3,7,11),]
pts_in_left=pip3d(verts_left,faces,grid_full)
in_left=which(pts_in_left==1)
grid_left=grid_full[in_left,]
# right cone
verts_right=pos_r[c(1,2,3,7,11),]
pts_in_right=pip3d(verts_right,faces,grid_left)
in_right=which(pts_in_right==1)
grid_both=grid_left[in_right,]
# points in joint cones
scatter3D(grid_both[,1],grid_both[,2],grid_both[,3], col="black",pch = 19, cex = 0.5,add = TRUE)

# now to figure out the "change in volume" function
range_bins=seq(0.5,13.5,by=1)# these are the boundaries - only going out to 13 m as there are edge effects
grid_pt_ranges=sqrt(grid_both[,1]^2+grid_both[,2]^2+grid_both[,3]^2)# euclidean distance to origin (left camera)
vol=vector(mode="numeric",length=13)
for (i in 1:13){
  # finsd points in range interval
  ind=which(grid_pt_ranges>range_bins[i] & grid_pt_ranges<=range_bins[i+1])
  # points that qualify times volume per point
  vol[i]=length(ind)*v_pt
}
cam_range=1:13
plot(cam_range,vol)

# fit 2 order polynomial to this
p_vol= coef(lm(vol ~cam_range +I(cam_range^2)))
y=p_vol[1]+p_vol[2]*cam_range+p_vol[3]*cam_range^2
lines(cam_range,y,col="green")

# now we create a detection function (logistic) to visualize probability
d50=8
sl=-2
d_range=seq(0,15,length.out = 100)
p_detect=1./(1+9^((d50-d_range)/sl))
plot(d_range,p_detect, type = "l")

# we compute the effective area based on volume and detection functions
params=c(p_vol,d50,sl)
# define the function
eff_vol_func<- function(x,params){(params[1]+params[2]*x+params[3]*x^2)*(1./(1+9^((params[4]-x)/params[5])))}
# for good measure, plot this
plot(d_range,eff_vol_func(d_range,params), type = "l")
# integrate
eff_vol=integrate(eff_vol_func,params,lower =0, upper =15)


# now to run the simulation
n = 1000 # number of fish
reps=1000
estimated_density=numeric(length=reps)
for (i in 1:reps){
  # now we add some random fish
  xf=runif(n)*25-12.5
  yf=runif(n)*15
  zf=runif(n)*25-12.5

  # now get the ones we can see
  # left cone
  fish_full=matrix(c(xf,yf,zf),length(zf),3)
  pts_in_left=pip3d(verts_left,faces,fish_full)
  in_left=which(pts_in_left==1)
  fish_left=matrix(c(xf[in_left],yf[in_left],zf[in_left]),length(xf[in_left]),3)
  # right cone
  pts_in_right=pip3d(verts_right,faces,fish_left)
  in_right=which(pts_in_right==1)
  xfb=fish_left[in_right,1]
  yfb=fish_left[in_right,2]
  zfb=fish_left[in_right,3]
  
  # estimate probability of detection
  fish_range=sqrt(xfb^2+yfb^2+zfb^2)
  prob=1./(1+9^((d50-fish_range)/sl))
  detected=numeric(length(prob))
  for (j in 1:length(prob)){detected[j]=rbinom(1,1,prob[j])}
  
  estimated_density[i]=sum(detected)/eff_vol$value
}


known_density = n/(25*25*15)# this is in fish/m3
hist(estimated_density)
abline(v=known_density,col="red")
abline(v=mean(estimated_density),col="green")
# take home - we are "close-ish" with this method 
# the change-in-volume function can be assumed to be "fixed" as it is dependent 
# on physical properties of the camera system.  The detection function, and by 
# extension the "effective volume", as it pertains to a specific species/imaging 
# condition is what is needed to be estimated from existing data 


