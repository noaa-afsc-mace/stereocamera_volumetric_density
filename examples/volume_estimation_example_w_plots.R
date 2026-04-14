# this example follows the volume estimation process step by step,including 3D plots for ilustrative purposes

rm(list = ls())

library(StereoCamVolume)
library(plotly)
library(patchwork)
# read in a calibration
Cal=read_calibration('accessory_functions/example_calibrations/matlab_caltech_example2.mat','matlab_caltech')

max_extent=4 # maximum we want to plot the field
scaling=1000# scaling value is 1000, because calibration is in units of mm
normT=max_extent*scaling # extent of visual field in mm, 4 m extent scaled by scaling :-)

# the following code block is adapted from the camera calibration toolbox
#http://robots.stanford.edu/cs223b04/JeanYvesCalib/

BASE_left = normT*(t(matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1),3,3)) %*%
                     t(matrix(c(1/Cal$Camera1$Intrinsic$focal_length$h, 0, 0, 0, 1/Cal$Camera1$Intrinsic$focal_length$v, 0, 0, 0, 1),3,3)) %*%
                     t(matrix(c(1, 0, -Cal$Camera1$Intrinsic$principal_point$h, 0, 1, -Cal$Camera1$Intrinsic$principal_point$v, 0, 0, 1),3,3)) %*%
                     t(matrix(c(0, Cal$Camera1$Intrinsic$image_size$width-1, Cal$Camera1$Intrinsic$image_size$width-1, 0, 0, 0, 0, Cal$Camera1$Intrinsic$image_size$height-1, Cal$Camera1$Intrinsic$image_size$height-1, 0, 1, 1, 1, 1, 1),5,3)))
IP_left  = matrix(rbind(BASE_left,matrix(0,3,5),BASE_left),3,15)

RmatLeft=t(matrix(c(Cal$Camera1$Extrinsic$R$r1c1,Cal$Camera1$Extrinsic$R$r2c1,Cal$Camera1$Extrinsic$R$r3c1,Cal$Camera1$Extrinsic$R$r1c2,Cal$Camera1$Extrinsic$R$r2c2,Cal$Camera1$Extrinsic$R$r3c2,Cal$Camera1$Extrinsic$R$r1c3,Cal$Camera1$Extrinsic$R$r2c3,Cal$Camera1$Extrinsic$R$r3c3),3,3))
TmatLeft=matrix(c(Cal$Camera1$Extrinsic$T$t1, Cal$Camera1$Extrinsic$T$t2, Cal$Camera1$Extrinsic$T$t3),3,1)

IP_left = RmatLeft %*% (IP_left - matrix(TmatLeft,3, 15))

BASE_right = normT*(t(matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1),3,3)) %*%
                      t(matrix(c(1/Cal$Camera2$Intrinsic$focal_length$h, 0, 0, 0, 1/Cal$Camera2$Intrinsic$focal_length$v, 0, 0, 0, 1),3,3)) %*%
                      t(matrix(c(1, 0, -Cal$Camera2$Intrinsic$principal_point$h, 0, 1, -Cal$Camera2$Intrinsic$principal_point$v, 0, 0, 1),3,3)) %*%
                      t(matrix(c(0, Cal$Camera2$Intrinsic$image_size$width-1, Cal$Camera2$Intrinsic$image_size$width-1, 0, 0, 0, 0, Cal$Camera2$Intrinsic$image_size$height-1, Cal$Camera2$Intrinsic$image_size$height-1, 0, 1, 1, 1, 1, 1),5,3)))
IP_right  = matrix(rbind(BASE_right,matrix(0,3,5),BASE_right),3,15)

RmatRight=t(matrix(c(Cal$Camera2$Extrinsic$R$r1c1,Cal$Camera2$Extrinsic$R$r2c1,Cal$Camera2$Extrinsic$R$r3c1,Cal$Camera2$Extrinsic$R$r1c2,Cal$Camera2$Extrinsic$R$r2c2,Cal$Camera2$Extrinsic$R$r3c2,Cal$Camera2$Extrinsic$R$r1c3,Cal$Camera2$Extrinsic$R$r2c3,Cal$Camera2$Extrinsic$R$r3c3),3,3))
TmatRight=matrix(c(Cal$Camera2$Extrinsic$T$t1, Cal$Camera2$Extrinsic$T$t2, Cal$Camera2$Extrinsic$T$t3),3,1)

IP_right = RmatRight %*% (IP_right - matrix(TmatRight,3, 15))

# some of this is redundant, which causes issues later on, so I'm paring it down here
IP_left=IP_left[,c(1,2,4,4,5,7,7,8,10,10,11,13,1,4,7,10,1)]
IP_right=IP_right[,c(1,2,4,4,5,7,7,8,10,10,11,13,1,4,7,10,1)]

#originally this is for a down facing camera, so we flip the z and y axes, and reverse the direction of the z axis.  Also, change from mm to m
xl=IP_left[1,]/scaling
yl=IP_left[3,]/scaling
zl=-IP_left[2,]/scaling

xr=IP_right[1,]/scaling
yr=IP_right[3,]/scaling
zr=-IP_right[2,]/scaling

# generate grid points
grid_size=0.2 # this is in m, smaller values = denser grid - more precision and more computation load.
# For this plotting we kep it large.  For computing volume, we would reduce it by a factor of 10
# making sure we envelop the entire view cones with points
xg_vec=seq(floor(min(c(xl,xr))*10)/10,ceiling(max(c(xl,xr))*10)/10,by=grid_size)
yg_vec=seq(floor(min(c(yl,yr))*10)/10,ceiling(max(c(yl,yr))*10)/10,by=grid_size)
zg_vec=seq(floor(min(c(zl,zr))*10)/10,ceiling(max(c(zl,zr))*10)/10,by=grid_size)
grid_frame=expand.grid(xg_vec,yg_vec,zg_vec)
colnames(grid_frame) <- c("x","y","z")
# remove 0 and less from y (z) axis
ind=which(grid_frame$y>0)
grid_frame=grid_frame[ind,]

grid_full=matrix(c(grid_frame$x,grid_frame$y,grid_frame$z),length(grid_frame$z),3)
# revert to original down facing view cone. Here the y axis is up-down.
grid_nat=matrix(c(grid_frame$x,-grid_frame$z,grid_frame$y),length(grid_frame$z),3)
grid_both=find_targets_in_view(Cal,scaling, grid_nat)
# now revert us back to the natural axis orientation
grid_both=matrix(c(grid_both[,1],grid_both[,3],-grid_both[,2]),length(grid_both[,1]),3)
# here we want to rotate grid and remove points that are sub surface
floor_position=c(-10,10,0.5)# this camera is tilted down 10 degrees, to the right ten degrees and located at a height of 0.5 m for the floor
a=floor_position[1]*pi/180 # tilt in radians
b=floor_position[2]*pi/180 # roll in radians
g=0 #  yaw not implemented
# assemble rotation matrix
Rot=matrix(data=c(cos(a)*cos(g),
                  cos(g)*sin(a)*sin(b)-sin(g)*cos(a),
                  cos(a)*sin(b)*cos(g)+sin(a)*sin(g),
                  sin(g)*cos(b),
                  sin(a)*sin(b)*sin(g)+cos(a)*cos(g),
                  sin(a)*sin(b)*cos(g)-cos(g)*sin(a),
                  -sin(b),
                  cos(b)*sin(a),
                  cos(a)*cos(b)),nrow=3,ncol=3,byrow=TRUE)
Trans=matrix(data=c(0,0,floor_position[3]),3,nrow(grid_both))# translation matrix for camera position
grid_both_rot=t(Rot%*%t(grid_both)+Trans)# rotate point grid
grid_both_rot=grid_both_rot[which(grid_both_rot[,3]>0),]# keep points above sea floor (e.g. positive points)

# lets rotate and translate the view cones too
Trans=matrix(data=c(0,0,floor_position[3]),3,length(xl))
left_field=t(Rot%*%t(matrix(c(xl,yl,zl),length(xl),3))+Trans)
right_field=t(Rot%*%t(matrix(c(xr,yr,zr),length(xr),3))+Trans)

#plot - this is a rotatable 3D figure of the view fields

# 1. Base Grid (Blue Markers)
fig <- plot_ly(x = ~grid_both_rot[,1], y = ~grid_both_rot[,2], z = ~grid_both_rot[,3],
               type = 'scatter3d',
               mode = 'markers',
               name = 'Grid Points',
               marker = list(size = 2, color = "blue", opacity = 0.5))

# 2. Right Camera (Green Line)
fig <- fig %>% add_trace(
  x = ~right_field[,1], # Standardized to first index for closing loops
  y = ~right_field[,2],
  z = ~right_field[,3],
  type = 'scatter3d',
  mode = 'lines',
  name = 'Right view cone',
  line = list(width = 5, color = "green")
)

# 3. Left Camera (Red Line)
fig <- fig %>% add_trace(
  x = ~left_field[,1],
  y = ~left_field[,2],
  z = ~left_field[,3],
  type = 'scatter3d',
  mode = 'lines',
  name = 'Left view cone',
  line = list(width = 5, color = "red")
)

# 4. Layout
fig <- fig %>% layout(
  scene = list(
    xaxis = list(title = "X Axis"),
    yaxis = list(title = "Y Axis"),
    zaxis = list(title = "Z Axis")
  ),
  showlegend = TRUE
)

#show(fig)
# another way to view the field
# dataPth=data.frame(cbind(x=grid_both_rot[,1],y=grid_both_rot[,2],z=grid_both_rot[,3]))
# p1=ggplot(dataPth, aes(x=x, y=y)) +
#   geom_point(shape = 16, size = 3,color="red") +
#   labs(title = "top", x = "X (m)", y = "Y (m)")  +
#   theme_bw()
# p2=ggplot(dataPth, aes(x=y, y=z)) +
#   geom_point(shape = 16, size = 3,color="red") +
#   labs(title = "side", x = "Y (m)", y = "Z (m)")  +
#   theme_bw()
# p3=ggplot(dataPth, aes(x=x, y=z)) +
#   geom_point(shape = 16, size = 3,color="red") +
#   labs(title = "front", x = "X (m)", y = "Z (m)")  +
#   theme_bw()
#
# (p1|p2|p3)


# now to figure out the "change in volume" function

camera_pos=c(0,0,floor_position[3])
grid_pt_ranges=sqrt((grid_both_rot[,1]-camera_pos[1])^2+
                      (grid_both_rot[,2]-camera_pos[2])^2+
                      (grid_both_rot[,3]-camera_pos[3])^2)# euclidean distance to origin (left camera)
range_bins=seq(min(grid_pt_ranges),round(max_extent)-0.5,length.out=25)# these are the boundaries - only going out to max extent -0.5 m as there are edge effects
range_bin_halfpoint=(range_bins[2]-range_bins[1])/2
c_vol=vector(mode="numeric",length=length(range_bins))
v_pt=grid_size^3
# The key here is that the bin width is equivalent to the unit, (e.g. 1 m).
# This gives us the change function.
for (i in 1:length(range_bins)){
  # find points in range interval
  ind=which(grid_pt_ranges<=range_bins[i]+range_bin_halfpoint)
  # number of points times volume per point
  c_vol[i]=length(ind)*v_pt
}

x0=min(grid_pt_ranges)# nearest range point
x=range_bins
# fit 3rd degree cumulative volume function (F) forced though nearest ranged point
pc_vol= as.numeric(coef(lm(c_vol ~-1 +x +I((x-x0)^2) +I((x-x0)^3))))
#make a plot
x=seq(0,max(range_bins),length.out=50)
y=pc_vol[1]*(x-x0)+pc_vol[2]*(x-x0)^2+pc_vol[3]*(x-x0)^3
plotdf=data.frame(range_bins=range_bins,c_vol=c_vol)
plotdf2=data.frame(cbind(x=x, y=y))
p1=ggplot() +
  geom_point(data=plotdf,aes(x=range_bins, y=c_vol, color = 'Point cloud est.'), shape = 3, size = 3) +
  geom_line(data=plotdf2, mapping=aes(x=x, y = y, color = 'Cum. vol. function'), linewidth = 0.5) +
  scale_color_manual(values = c('black', 'black')) +
  guides(color = guide_legend("")) +
  # Add the arrow
  annotate("segment", x = x0, y = 2, xend = x0, yend = 0,
           arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  # Add accompanying text
  annotate("text", x = x0, y = 2.5, label = expression("x"[0]), color = "black") +
  labs(title = "", x = "range from camera (m)", y = expression("cumulative volume (m"^"3)"),color="Legend")  +
  theme_bw()+ theme(legend.position = "top")
p1


#
vol_func<- function(x){ifelse(x<=x0,0,2*pc_vol[2]*(x-x0)+3*pc_vol[3]*(x-x0)^2)}
# fit 2 order polynomial to this

y=vol_func(x)
ind1=which(x==1.5)
ind2=which(x==2.5)
ind=c(seq(1,ind1),seq(ind1,ind2),seq(ind2,length(x)))
ya=y[ind]

ya[seq(1,ind1)]=0
ya[seq(ind2+2,length(ya))]=0
x=x[ind]
y=vol_func(x)
plotdf2=data.frame(cbind(x=x, y=y))
p2=ggplot(plotdf2, aes(x=x, y=y)) +
  geom_line(aes(y = y, color = "Vol. change function"), linewidth = 0.5) +
  geom_area(aes(y = ya, fill = "Volume"), alpha = 0.4) +

  scale_color_manual(name = "", values = c("Vol. change function" = "black")) +
  scale_fill_manual(name = "", values = c("Volume" = "gray45")) +
  labs(title = "", x = "range from camera (m)", y = expression("change in volume (m"^"3)"),color="Legend")  +

  theme_bw()+ theme(legend.position = "top")
# show(p2)


(p1|p2)




