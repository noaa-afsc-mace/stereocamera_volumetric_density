rm(list = ls())
library(yaml)
library(ptinpoly)
library(plot3D)
cal='calibrations/example_calibration_fromopencv.yml'
Cal=yaml.load(read_yaml(cal))
max_extent=2
grid_size=0.5
include_distortion=FALSE
  # code below is adopted from the camera calibration toolbox
  #http://robots.stanford.edu/cs223b04/JeanYvesCalib/
  # cited as Bouguet, J.Y., 2008. Camera calibration toolbox for Matlab [online]. [Available from http://vision.caltech.edu/bouguetj/calib_doc/index.html (accessed September 2008)].

  # this code sets up the view "cones"
  normT=max_extent*1000 # extent of visual field in mm

  BASE_left = normT*(t(matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1),3,3)) %*%
                       t(matrix(c(1/Cal$Camera1$Intrinsic$focal_point$h, 0, 0, 0, 1/Cal$Camera1$Intrinsic$focal_point$v, 0, 0, 0, 1),3,3)) %*%
                       t(matrix(c(1, 0, -Cal$Camera1$Intrinsic$principal_point$h, 0, 1, -Cal$Camera1$Intrinsic$principal_point$v, 0, 0, 1),3,3)) %*%
                       t(matrix(c(0, Cal$Camera1$Intrinsic$image_size$width-1, Cal$Camera1$Intrinsic$image_size$width-1, 0, 0, 0, 0, Cal$Camera1$Intrinsic$image_size$height-1, Cal$Camera1$Intrinsic$image_size$height-1, 0, 1, 1, 1, 1, 1),5,3)))
  IP_left  = matrix(rbind(BASE_left,matrix(0,3,5),BASE_left),3,15)

  RmatLeft=t(matrix(c(Cal$Camera1$Extrinsic$R$r1c1,Cal$Camera1$Extrinsic$R$r2c1,Cal$Camera1$Extrinsic$R$r3c1,Cal$Camera1$Extrinsic$R$r1c2,Cal$Camera1$Extrinsic$R$r2c2,Cal$Camera1$Extrinsic$R$r3c2,Cal$Camera1$Extrinsic$R$r1c3,Cal$Camera1$Extrinsic$R$r2c3,Cal$Camera1$Extrinsic$R$r3c3),3,3))
  TmatLeft=matrix(c(Cal$Camera1$Extrinsic$T$t1, Cal$Camera1$Extrinsic$T$t2, Cal$Camera1$Extrinsic$T$t3),3,1)

  IP_left = RmatLeft %*% (IP_left - matrix(TmatLeft,3, 15))

  BASE_right = normT*(t(matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1),3,3)) %*%
                       t(matrix(c(1/Cal$Camera2$Intrinsic$focal_point$h, 0, 0, 0, 1/Cal$Camera2$Intrinsic$focal_point$v, 0, 0, 0, 1),3,3)) %*%
                       t(matrix(c(1, 0, -Cal$Camera2$Intrinsic$principal_point$h, 0, 1, -Cal$Camera2$Intrinsic$principal_point$v, 0, 0, 1),3,3)) %*%
                       t(matrix(c(0, Cal$Camera2$Intrinsic$image_size$width-1, Cal$Camera2$Intrinsic$image_size$width-1, 0, 0, 0, 0, Cal$Camera2$Intrinsic$image_size$height-1, Cal$Camera2$Intrinsic$image_size$height-1, 0, 1, 1, 1, 1, 1),5,3)))
  IP_right  = matrix(rbind(BASE_right,matrix(0,3,5),BASE_right),3,15)

  RmatRight=t(matrix(c(Cal$Camera2$Extrinsic$R$r1c1,Cal$Camera2$Extrinsic$R$r2c1,Cal$Camera2$Extrinsic$R$r3c1,Cal$Camera2$Extrinsic$R$r1c2,Cal$Camera2$Extrinsic$R$r2c2,Cal$Camera2$Extrinsic$R$r3c2,Cal$Camera2$Extrinsic$R$r1c3,Cal$Camera2$Extrinsic$R$r2c3,Cal$Camera2$Extrinsic$R$r3c3),3,3))
  TmatRight=matrix(c(Cal$Camera2$Extrinsic$T$t1, Cal$Camera2$Extrinsic$T$t2, Cal$Camera2$Extrinsic$T$t3),3,1)

  IP_right = RmatRight %*% (IP_right - matrix(TmatRight,3, 15))

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

  lines3D(xl,yl,zl,xlab = "x (m)", ylab = "y (m)", zlab = "z (m)", col="red",ticktype = "detailed")
  lines3D(xr,yr,zr, col="green",add = TRUE)
  # generate grid points
  # making sure we envelop the entire view cones with points
  xg_vec=seq(floor(min(c(xl,xr))*10)/10,ceiling(max(c(xl,xr))*10)/10,by=grid_size)
  yg_vec=seq(floor(min(c(yl,yr))*10)/10,ceiling(max(c(yl,yr))*10)/10,by=grid_size)
  zg_vec=seq(floor(min(c(zl,zr))*10)/10,ceiling(max(c(zl,zr))*10)/10,by=grid_size)
  grid_frame=expand.grid(xg_vec,yg_vec,zg_vec)
  colnames(grid_frame) <- c("x","y","z")
  # remove 0 and less from y (z) axis
  ind=which(grid_frame$y>0)
  grid_frame=grid_frame[ind,]
  # find ones that are within both cones
  pos_l=matrix(c(xl,yl,zl),length(xl),3)
  pos_r=cbind(xr,yr,zr)



  grid_full=matrix(c(grid_frame$x,grid_frame$y,grid_frame$z),length(grid_frame$z),3)
  grid_nat=matrix(c(grid_frame$x,-grid_frame$z,grid_frame$y),length(grid_frame$z),3)
  # pixel projection method
  #adopted from https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html
  #left projection (no translation/rotation)
  grid_nat_rot=t(RmatLeft %*% t(grid_nat)+matrix(TmatLeft/1000,3, nrow(grid_nat)))
  xp=grid_nat_rot[,1]/grid_nat_rot[,3]
  yp=grid_nat_rot[,2]/grid_nat_rot[,3]
  r=sqrt(xp^2+yp^2)
  xf=xp*(1+Cal$Camera1$Intrinsic$radial_distortion$k1*r^2+Cal$Camera1$Intrinsic$radial_distortion$k2*r^4+Cal$Camera1$Intrinsic$radial_distortion$k3*r^6)+2*xp*yp*+Cal$Camera1$Intrinsic$tangential_distortion$p1+Cal$Camera1$Intrinsic$tangential_distortion$p2*(r^2+2*xp^2)
  yf=yp*(1+Cal$Camera1$Intrinsic$radial_distortion$k1*r^2+Cal$Camera1$Intrinsic$radial_distortion$k2*r^4+Cal$Camera1$Intrinsic$radial_distortion$k3*r^6)+2*xp*yp*+Cal$Camera1$Intrinsic$tangential_distortion$p2+Cal$Camera1$Intrinsic$tangential_distortion$p1*(r^2+2*yp^2)
  u=Cal$Camera1$Intrinsic$focal_point$h*xf+Cal$Camera1$Intrinsic$principal_point$h
  v=Cal$Camera1$Intrinsic$focal_point$v*yf+Cal$Camera1$Intrinsic$principal_point$v
  in_left=which(u>0 & u<Cal$Camera1$Intrinsic$image_size$width & v>0 & v<Cal$Camera1$Intrinsic$image_size$height)


  points3D(grid_nat[in_left,1],grid_nat[in_left,3],-grid_nat[in_left,2],col="blue")

  grid_left=grid_nat[in_left,]
  #right projection (impose translation/rotation)
  grid_left_rot=t(RmatRight %*% t(grid_left)+matrix(TmatRight/1000,3, nrow(grid_left)))

  xp=grid_left_rot[,1]/grid_left_rot[,3]
  yp=grid_left_rot[,2]/grid_left_rot[,3]
  r=sqrt(xp^2+yp^2)
  xf=xp*(1+Cal$Camera2$Intrinsic$radial_distortion$k1*r^2+Cal$Camera2$Intrinsic$radial_distortion$k2*r^4+Cal$Camera2$Intrinsic$radial_distortion$k3*r^6)+2*xp*yp*+Cal$Camera2$Intrinsic$tangential_distortion$p1+Cal$Camera2$Intrinsic$tangential_distortion$p2*(r^2+2*xp^2)
  yf=yp*(1+Cal$Camera2$Intrinsic$radial_distortion$k1*r^2+Cal$Camera2$Intrinsic$radial_distortion$k2*r^4+Cal$Camera2$Intrinsic$radial_distortion$k3*r^6)+2*xp*yp*+Cal$Camera2$Intrinsic$tangential_distortion$p2+Cal$Camera2$Intrinsic$tangential_distortion$p1*(r^2+2*yp^2)
  u=Cal$Camera2$Intrinsic$focal_point$h*xf+Cal$Camera2$Intrinsic$principal_point$h
  v=Cal$Camera2$Intrinsic$focal_point$v*yf+Cal$Camera2$Intrinsic$principal_point$v
  in_right=which(u>0 & u<Cal$Camera2$Intrinsic$image_size$width & v>0 & v<Cal$Camera2$Intrinsic$image_size$height)
  grid_both=grid_left[in_right,]
  points3D(grid_both[,1],grid_both[,3],-grid_both[,2],col="red")



  # faces=t(matrix(c(1,2,3,1,4,3,4,5,3,5,2,3,1,2,4,4,5,2),3,6))
  # verts_left=pos_l[c(1,9,2,3,6),]
  # pts_in_left=ptinpoly::pip3d(verts_left,faces,grid_full)
  # in_left=which(pts_in_left==1)
  # grid_left=grid_full[in_left,]
  # # right cone
  # verts_right=pos_r[c(1,9,2,3,6),]
  # pts_in_right=ptinpoly::pip3d(verts_right,faces,grid_left)
  # in_right=which(pts_in_right==1)
  # grid_both=grid_left[in_right,]

  # now to figure out the "change in volume" function
  range_bins=seq(0.5,7.5,by=1)# these are the boundaries - only going out to 13 m as there are edge effects
  grid_pt_ranges=sqrt(grid_both[,1]^2+grid_both[,2]^2+grid_both[,3]^2)# euclidean distance to origin (left camera)
  vol=vector(mode="numeric",length=7)
  v_pt=grid_size^3
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
  p_vol= as.numeric(coef(lm(vol ~cam_range +I(cam_range^2))))
  disp_range=seq(0,7,length.out = 100)
  y=p_vol[1]+p_vol[2]*disp_range+p_vol[3]*disp_range^2

  # create the function object for itegration
  vol_func<- function(x){p_vol[1]+p_vol[2]*x+p_vol[3]*x^2}


