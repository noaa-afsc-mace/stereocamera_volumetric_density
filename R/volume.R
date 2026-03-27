#' Get volume function (change in joint imaging volume with range) from a calibration.
#' @param cal A yaml file (full path) from which to read the calibration values.
#' @param max_extent distance in m determining range extent of volume field.
#' @param grid_size grid spacing in m determining the density of the point cloud.
#' @param plotting boolean flag to indicate whether to make plots
#' @param units this parameter sets the general volumetric units of the analysis (cubic meters ("m3") for larger systems, or liter ("l") for smaller ones)
#' @param floor_position a set of three parameters indicating two angles (roll, then tilt) of the seafloor plane and theheight of the camera from that plane
#' @return A list containing the function vol_func and a vector of coefficients p_vol
#' @references Bouget (2008) Camera calibration toolbox for Matlab. Available from http://vision.caltech.edu/bouguetj/calib_doc/index.html (accessed September 2008)].
#' @details to do
#' @export
get_vol_func <- function(Cal, max_extent=8, grid_size=0.1, plotting=FALSE, units='m3', floor_position=c(NA,NA,NA)){
  # code below is adapted from the camera calibration toolbox
  #http://robots.stanford.edu/cs223b04/JeanYvesCalib/
  # cited as Bouguet, J.Y., 2008. Camera calibration toolbox for Matlab [online]. [Available from http://vision.caltech.edu/bouguetj/calib_doc/index.html (accessed September 2008)].
  # Calibration (Cal) should be a list of lists from reading in standard yml file
  # this code sets up the view "cones"
  if (units=='m3'){
    scaling=1000
    xlab='range from camera (m)'

  }
  else if (units=='l'){
    scaling=100
    xlab='range from camera (dm)'
  }
  normT=max_extent*scaling # extent of visual field in mm
  # rough estimation of viewfields for point cloud generation
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
  IP_left=IP_left[,c(1,2,4,4,5,7,7,8,10,10,11,13,1,4,7,10)]
  IP_right=IP_right[,c(1,2,4,4,5,7,7,8,10,10,11,13,1,4,7,10)]

  #originally this is for a down facing camera, so we flip the z and y axes, and reverse the direction of the z axis.  Also, change from mm to m
  xl=IP_left[1,]/scaling
  yl=IP_left[3,]/scaling
  zl=-IP_left[2,]/scaling

  xr=IP_right[1,]/scaling
  yr=IP_right[3,]/scaling
  zr=-IP_right[2,]/scaling


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
  # revret to original down facing view cone. Here the y axis is up-down.
  grid_nat=matrix(c(grid_frame$x,-grid_frame$z,grid_frame$y),length(grid_frame$z),3)
  grid_both=find_targets_in_view(Cal,scaling, grid_nat)
  # now revert us back to the natural axis orientation
  grid_both=matrix(c(grid_both[,1],grid_both[,3],-grid_both[,2]),length(grid_both[,1]),3)

  if (!is.na(floor_position[1])) {
    # here we want to rotate grid and remove points that are sub surface
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
    Trans=matrix(data=c(0,0,floor_position[3]),3,nrow(grid_both))
    grid_both_rot=t(Rot%*%t(grid_both)+Trans)# rotate point grid
    grid_both_rot=grid_both_rot[which(grid_both_rot[,3]>0),]# keep points above sea floor (e.g. positive points)
    Rotinv=solve(Rot)
    Trans=matrix(data=c(0,0,floor_position[3]),3,nrow(grid_both_rot))
    grid_both=t(Rotinv%*%t(grid_both_rot)-Trans)
  }


  # now to figure out the "change in volume" function
  range_bins=seq(0.5,round(max_extent)-0.5,by=1)# these are the boundaries - only going out to max_extent m as there are edge effects
  grid_pt_ranges=sqrt(grid_both[,1]^2+grid_both[,2]^2+grid_both[,3]^2)# euclidean distance to origin (left camera)
  vol=vector(mode="numeric",length=round(max_extent)-1)
  v_pt=grid_size^3
  # The key here is that the bin width is equivalent to the unit, (e.g. 1 m).
  # This gives us the change function.
  for (i in 1:round(max_extent)-1){
    # find points in range interval
    ind=which(grid_pt_ranges>range_bins[i] & grid_pt_ranges<=range_bins[i+1])
    # number of points times volume per point
    vol[i]=length(ind)*v_pt
  }
  range_centers=0:(round(max_extent)-1)
  vol=c(0,vol)
  # fit 2 order polynomial to this
  p_vol= as.numeric(coef(lm(vol ~range_centers +I(range_centers^2))))

  # create the function object for itegration
  vol_func<- function(x){p_vol[1]+p_vol[2]*x+p_vol[3]*x^2}

  # plot if asked to
  if (plotting){
    plot(range_centers,vol,xlab=xlab,ylab='change in volume')
    x=seq(0, max(range_centers),length.out=100)
    y=vol_func(x)
    lines(x,y,col="red")
  }

  return(list(vol_func=vol_func, p_vol=p_vol, vol=vol, range_centers=range_centers))
}


#' @param cal A yaml file (full path) from which to read the calibration values.
#' @param target_positions 3d coordinates in appropriate units.
#' @param scaling this relates calibration units to real world.
#' @return target positions within both viewing cones
#' @references Williams, K., Rooper, C.N., De Robertis, A., Levine, M. and Towler, R., 2018. A method for computing volumetric fish density using stereo cameras. Journal of Experimental Marine Biology and Ecology, 508, pp.21-26.
#' @details to do
#' @export
find_targets_in_view <- function(Cal, scaling, target_positions){
  RmatLeft=t(matrix(c(Cal$Camera1$Extrinsic$R$r1c1,Cal$Camera1$Extrinsic$R$r2c1,Cal$Camera1$Extrinsic$R$r3c1,Cal$Camera1$Extrinsic$R$r1c2,Cal$Camera1$Extrinsic$R$r2c2,Cal$Camera1$Extrinsic$R$r3c2,Cal$Camera1$Extrinsic$R$r1c3,Cal$Camera1$Extrinsic$R$r2c3,Cal$Camera1$Extrinsic$R$r3c3),3,3))
  TmatLeft=matrix(c(Cal$Camera1$Extrinsic$T$t1, Cal$Camera1$Extrinsic$T$t2, Cal$Camera1$Extrinsic$T$t3),3,1)
  RmatRight=t(matrix(c(Cal$Camera2$Extrinsic$R$r1c1,Cal$Camera2$Extrinsic$R$r2c1,Cal$Camera2$Extrinsic$R$r3c1,Cal$Camera2$Extrinsic$R$r1c2,Cal$Camera2$Extrinsic$R$r2c2,Cal$Camera2$Extrinsic$R$r3c2,Cal$Camera2$Extrinsic$R$r1c3,Cal$Camera2$Extrinsic$R$r2c3,Cal$Camera2$Extrinsic$R$r3c3),3,3))
  TmatRight=matrix(c(Cal$Camera2$Extrinsic$T$t1, Cal$Camera2$Extrinsic$T$t2, Cal$Camera2$Extrinsic$T$t3),3,1)
  # pixel projection method
  #adopted from https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html
  #left projection (impose translation/rotation, although for matlab / opencv this doesn't change values
  target_positions=t(RmatLeft %*% t(target_positions)+matrix(TmatLeft/scaling,3, nrow(target_positions)))
  xp=target_positions[,1]/target_positions[,3]
  yp=target_positions[,2]/target_positions[,3]
  r=sqrt(xp^2+yp^2)
  xf=xp*(1+Cal$Camera1$Intrinsic$radial_distortion$k1*r^2+Cal$Camera1$Intrinsic$radial_distortion$k2*r^4+Cal$Camera1$Intrinsic$radial_distortion$k3*r^6)+2*xp*yp*+Cal$Camera1$Intrinsic$tangential_distortion$p1+Cal$Camera1$Intrinsic$tangential_distortion$p2*(r^2+2*xp^2)
  yf=yp*(1+Cal$Camera1$Intrinsic$radial_distortion$k1*r^2+Cal$Camera1$Intrinsic$radial_distortion$k2*r^4+Cal$Camera1$Intrinsic$radial_distortion$k3*r^6)+2*xp*yp*+Cal$Camera1$Intrinsic$tangential_distortion$p2+Cal$Camera1$Intrinsic$tangential_distortion$p1*(r^2+2*yp^2)
  u=Cal$Camera1$Intrinsic$focal_length$h*xf+Cal$Camera1$Intrinsic$principal_point$h
  v=Cal$Camera1$Intrinsic$focal_length$v*yf+Cal$Camera1$Intrinsic$principal_point$v
  in_left=which(u>0 & u<Cal$Camera1$Intrinsic$image_size$width & v>0 & v<Cal$Camera1$Intrinsic$image_size$height)
  target_positions_left=target_positions[in_left,]

  #right projection (impose translation/rotation)
  grid_left_rot=t(RmatRight %*% t(target_positions_left)+matrix(TmatRight/scaling,3, nrow(target_positions_left)))
  xp=grid_left_rot[,1]/grid_left_rot[,3]
  yp=grid_left_rot[,2]/grid_left_rot[,3]
  r=sqrt(xp^2+yp^2)
  xf=xp*(1+Cal$Camera2$Intrinsic$radial_distortion$k1*r^2+Cal$Camera2$Intrinsic$radial_distortion$k2*r^4+Cal$Camera2$Intrinsic$radial_distortion$k3*r^6)+2*xp*yp*+Cal$Camera2$Intrinsic$tangential_distortion$p1+Cal$Camera2$Intrinsic$tangential_distortion$p2*(r^2+2*xp^2)
  yf=yp*(1+Cal$Camera2$Intrinsic$radial_distortion$k1*r^2+Cal$Camera2$Intrinsic$radial_distortion$k2*r^4+Cal$Camera2$Intrinsic$radial_distortion$k3*r^6)+2*xp*yp*+Cal$Camera2$Intrinsic$tangential_distortion$p2+Cal$Camera2$Intrinsic$tangential_distortion$p1*(r^2+2*yp^2)
  u=Cal$Camera2$Intrinsic$focal_length$h*xf+Cal$Camera2$Intrinsic$principal_point$h
  v=Cal$Camera2$Intrinsic$focal_length$v*yf+Cal$Camera2$Intrinsic$principal_point$v
  in_right=which(u>0 & u<Cal$Camera2$Intrinsic$image_size$width & v>0 & v<Cal$Camera2$Intrinsic$image_size$height)
  target_positions_both=target_positions_left[in_right,]
  return(target_positions_both)

}


#' @param filename path to calibration file. For seagis files, it needs to be two files, one for left and one for right in that order
#' @param method this flag refers to which type of calibration file is being provided. Curent options are matlab_caltech, opencv, and seagis.
#' @return calibration structure
#' @references
#' @details to do
#' @export
read_calibration <- function(filename,method='matlab_caltech'){
  if (method=='matlab_caltech'){
    library(R.matlab)
    Calraw = readMat(filename)
    Camera1=list()
    Camera1$Extrinsic=list()
    # for matlab and opencv, the Camera1 camera coordinate system is considered the prime.  For SeaGIS, both Camera1 and Camera2 cameras are rotated from a coordinate system set at the center of the two cameras.
    Camera1$Extrinsic$R=list(r1c1=1,r1c2=0,r1c3=0,
                             r2c1=0,r2c2=1,r2c3=0,
                             r3c1=0,r3c2=0,r3c3=1)
    Camera1$Extrinsic$T=list(t1=0,t2=0,t3=0)

    Camera1$Intrinsic=list()
    Camera1$Intrinsic$focal_length=list(h=Calraw$fc.left[1],v=Calraw$fc.left[2])
    Camera1$Intrinsic$principal_point=list(h=Calraw$cc.left[1],v=Calraw$cc.left[2])
    Camera1$Intrinsic$image_size=list(width=Calraw$image.size.left[1],height=Calraw$image.size.left[2])
    Camera1$Intrinsic$radial_distortion=list(k1=Calraw$kc.left[1],k2=Calraw$kc.left[2],k3=Calraw$kc.left[5])
    Camera1$Intrinsic$tangential_distortion=list(p1=Calraw$kc.left[3],p2=Calraw$kc.left[4])

    Camera2=list()
    Camera2$Extrinsic=list()
    Camera2$Extrinsic$R=list(r1c1=Calraw$R[1,1],r1c2=Calraw$R[1,2],r1c3=Calraw$R[1,3],
                             r2c1=Calraw$R[2,1],r2c2=Calraw$R[2,2],r2c3=Calraw$R[2,3],
                             r3c1=Calraw$R[3,1],r3c2=Calraw$R[3,2],r3c3=Calraw$R[3,3])
    Camera2$Extrinsic$T=list(t1=Calraw$T[1],t2=Calraw$T[2],t3=Calraw$T[3])

    Camera2$Intrinsic=list()
    Camera2$Intrinsic$focal_length=list(h=Calraw$fc.right[1],v=Calraw$fc.right[2])
    Camera2$Intrinsic$principal_point=list(h=Calraw$cc.right[1],v=Calraw$cc.right[2])
    Camera2$Intrinsic$image_size=list(width=Calraw$image.size.right[1],height=Calraw$image.size.right[2])
    Camera2$Intrinsic$radial_distortion=list(k1=Calraw$kc.right[1],k2=Calraw$kc.right[2],k3=Calraw$kc.right[5])
    Camera2$Intrinsic$tangential_distortion=list(p1=Calraw$kc.right[3],p2=Calraw$kc.right[4])
    return(list(Camera1=Camera1, Camera2=Camera2))
  }

  else if (method=='seagis'){
    if (length(filename)!=2){
      stop("need both Camera2 and Camera1 files!")
    }
    left_file_name=filename[1]
    left_file_data=read.delim(toString(left_file_name),header = TRUE,skip = 4,sep = "\t")

    Camera1=list()
    Camera1$Extrinsic=list()
    # for matlab and opencv, the Camera1 camera coordinate system is considered the prime.  For SeaGIS, both Camera1 and Camera2 cameras are rotated from a coordinate system set at the center of the two cameras.
    # convert from rotation vector to matrix (Rodrigues)
    om=matrix(c(as.numeric(left_file_data$Data[which(left_file_data$Name=="Omega")]),
                as.numeric(left_file_data$Data[which(left_file_data$Name=="Phi")]),
                as.numeric(left_file_data$Data[which(left_file_data$Name=="Kappa")])),3,1)*pi/180*-1
    theta=norm(om,type = "2")
    omn=om/theta
    K = t(matrix(c(0, -omn[3,1], omn[2,1],omn[3,1], 0, -omn[1,1],-omn[2,1], omn[1,1], 0),3,3))
    R = diag(3) + sin(theta) * K + (1 - cos(theta)) * K %*% K

    Camera1$Extrinsic$R=list(r1c1=R[1,1],r1c2=R[1,2],r1c3=R[1,3],
                             r2c1=R[2,1],r2c2=R[2,2],r2c3=R[2,3],
                             r3c1=R[3,1],r3c2=R[3,2],r3c3=R[3,3])
    Camera1$Extrinsic$T=list(t1=as.numeric(left_file_data$Data[which(left_file_data$Name=="Camera X")]),
                             t2=as.numeric(left_file_data$Data[which(left_file_data$Name=="Camera Y")]),
                             t3=as.numeric(left_file_data$Data[which(left_file_data$Name=="Camera Z")]))

    Camera1$Intrinsic=list()
    fl_x=as.numeric(left_file_data$Data[which(left_file_data$Name=="Focal length")])/as.numeric(left_file_data$Data[which(left_file_data$Name=="X Pixel size")])
    fl_y=as.numeric(left_file_data$Data[which(left_file_data$Name=="Focal length")])/as.numeric(left_file_data$Data[which(left_file_data$Name=="Y Pixel size")])
    Camera1$Intrinsic$focal_length=list(h=fl_x,v=fl_y)
    pp_x=as.numeric(left_file_data$Data[which(left_file_data$Name=="Image columns")])/2+as.numeric(left_file_data$Data[which(left_file_data$Name=="X PP Offset")])/as.numeric(left_file_data$Data[which(left_file_data$Name=="X Pixel size")])
    pp_y=as.numeric(left_file_data$Data[which(left_file_data$Name=="Image rows")])/2+as.numeric(left_file_data$Data[which(left_file_data$Name=="Y PP Offset")])/as.numeric(left_file_data$Data[which(left_file_data$Name=="Y Pixel size")])
    Camera1$Intrinsic$principal_point=list(h=pp_x,v=pp_y)
    Camera1$Intrinsic$image_size=list(width=as.numeric(left_file_data$Data[which(left_file_data$Name=="Image columns")]),height=as.numeric(left_file_data$Data[which(left_file_data$Name=="Image rows")]))
    k1=as.numeric(left_file_data$Data[which(left_file_data$Name=="Radial distortion (K3)")])
    k2=as.numeric(left_file_data$Data[which(left_file_data$Name=="Radial distortion (K5)")])
    k3=as.numeric(left_file_data$Data[which(left_file_data$Name=="Radial distortion (K7)")])
    p1=as.numeric(left_file_data$Data[which(left_file_data$Name=="Decentring distortion (P1)")])
    p2=as.numeric(left_file_data$Data[which(left_file_data$Name=="Decentring distortion (P2)")])

    Camera1$Intrinsic$radial_distortion=list(k1=k1,k2=k2,k3=k3)
    Camera1$Intrinsic$tangential_distortion=list(p1=p1,p2=p2)

    right_file_name=filename[2]
    right_file_data=read.delim(toString(right_file_name),header = TRUE,skip = 4,sep = "\t")
    # convert from rotation vector to matrix (Rodrigues)
    om=matrix(c(as.numeric(right_file_data$Data[which(right_file_data$Name=="Omega")]),
                as.numeric(right_file_data$Data[which(right_file_data$Name=="Phi")]),
                as.numeric(right_file_data$Data[which(right_file_data$Name=="Kappa")])),3,1)*pi/180*-1
    theta=norm(om,type = "2")
    omn=om/theta
    K = t(matrix(c(0, -omn[3,1], omn[2,1],omn[3,1], 0, -omn[1,1],-omn[2,1], omn[1,1], 0),3,3))
    R = diag(3) + sin(theta) * K + (1 - cos(theta)) * K %*% K


    Camera2=list()
    Camera2$Extrinsic$R=list(r1c1=R[1,1],r1c2=R[1,2],r1c3=R[1,3],
                             r2c1=R[2,1],r2c2=R[2,2],r2c3=R[2,3],
                             r3c1=R[3,1],r3c2=R[3,2],r3c3=R[3,3])
    Camera2$Extrinsic$T=list(t1=as.numeric(right_file_data$Data[which(right_file_data$Name=="Camera X")]),
                             t2=as.numeric(right_file_data$Data[which(right_file_data$Name=="Camera Y")]),
                             t3=as.numeric(right_file_data$Data[which(right_file_data$Name=="Camera Z")]))

    Camera2$Intrinsic=list()
    fl_x=as.numeric(right_file_data$Data[which(right_file_data$Name=="Focal length")])/as.numeric(right_file_data$Data[which(right_file_data$Name=="X Pixel size")])
    fl_y=as.numeric(right_file_data$Data[which(right_file_data$Name=="Focal length")])/as.numeric(right_file_data$Data[which(right_file_data$Name=="Y Pixel size")])
    Camera2$Intrinsic$focal_length=list(h=fl_x,v=fl_y)
    pp_x=as.numeric(right_file_data$Data[which(right_file_data$Name=="Image columns")])/2+as.numeric(right_file_data$Data[which(right_file_data$Name=="X PP Offset")])/as.numeric(right_file_data$Data[which(right_file_data$Name=="X Pixel size")])
    pp_y=as.numeric(right_file_data$Data[which(right_file_data$Name=="Image rows")])/2+as.numeric(right_file_data$Data[which(right_file_data$Name=="Y PP Offset")])/as.numeric(right_file_data$Data[which(right_file_data$Name=="Y Pixel size")])
    Camera2$Intrinsic$principal_point=list(h=pp_x,v=pp_y)
    Camera2$Intrinsic$image_size=list(width=as.numeric(right_file_data$Data[which(right_file_data$Name=="Image columns")]),height=as.numeric(right_file_data$Data[which(right_file_data$Name=="Image rows")]))
    k1=as.numeric(right_file_data$Data[which(right_file_data$Name=="Radial distortion (K3)")])
    k2=as.numeric(right_file_data$Data[which(right_file_data$Name=="Radial distortion (K5)")])
    k3=as.numeric(right_file_data$Data[which(right_file_data$Name=="Radial distortion (K7)")])
    p1=as.numeric(right_file_data$Data[which(right_file_data$Name=="Decentring distortion (P1)")])
    p2=as.numeric(right_file_data$Data[which(right_file_data$Name=="Decentring distortion (P2)")])

    Camera2$Intrinsic$radial_distortion=list(k1=k1,k2=k2,k3=k3)
    Camera2$Intrinsic$tangential_distortion=list(p1=p1,p2=p2)

    return(list(Camera2=Camera2, Camera1=Camera1))
  }
  else if (method=='opencv'){
    library(reticulate)
    np <- import("numpy")
    Calraw=np$load(filename)
    C1Matrix=Calraw['cameraMatrixL']
    C2Matrix=Calraw['cameraMatrixR']
    C1distCoeffs=Calraw['distCoeffsL']
    C2distCoeffs=Calraw['distCoeffsR']
    C2R=Calraw['R']
    C2T=Calraw['T']
    C1imgSize=Calraw['imageSizeL']
    C2imgSize=Calraw['imageSizeR']

    Camera1=list()
    Camera1$Extrinsic=list()
    # for matlab and opencv, the Camera1 camera coordinate system is considered the prime.  For SeaGIS, both Camera1 and Camera2 cameras are rotated from a coordinate system set at the center of the two cameras.
    Camera1$Extrinsic$R=list(r1c1=1,r1c2=0,r1c3=0,
                             r2c1=0,r2c2=1,r2c3=0,
                             r3c1=0,r3c2=0,r3c3=1)
    Camera1$Extrinsic$T=list(t1=0,t2=0,t3=0)

    Camera1$Intrinsic=list()
    Camera1$Intrinsic$focal_length=list(h=C1Matrix[1,1],v=C1Matrix[2,2])
    Camera1$Intrinsic$principal_point=list(h=C1Matrix[1,3],v=C1Matrix[2,3])
    Camera1$Intrinsic$image_size=list(width=C1imgSize[1],height=C1imgSize[2])
    Camera1$Intrinsic$radial_distortion=list(k1=C1distCoeffs[1,1],k2=C1distCoeffs[1,2],k3=C1distCoeffs[1,5])
    Camera1$Intrinsic$tangential_distortion=list(p1=C1distCoeffs[1,3],p2=C1distCoeffs[1,4])

    Camera2=list()
    Camera2$Extrinsic=list()
    Camera2$Extrinsic$R=list(r1c1=C2R[1,1],r1c2=C2R[1,2],r1c3=C2R[1,3],
                             r2c1=C2R[2,1],r2c2=C2R[2,2],r2c3=C2R[2,3],
                             r3c1=C2R[3,1],r3c2=C2R[3,2],r3c3=C2R[3,3])
    Camera2$Extrinsic$T=list(t1=C2T[1],t2=C2T[2],t3=C2T[3])

    Camera2$Intrinsic=list()
    Camera2$Intrinsic$focal_length=list(h=C2Matrix[1,1],v=C2Matrix[2,2])
    Camera2$Intrinsic$principal_point=list(h=C2Matrix[1,3],v=C2Matrix[2,3])
    Camera2$Intrinsic$image_size=list(width=C2imgSize[1],height=C2imgSize[2])
    Camera2$Intrinsic$radial_distortion=list(k1=C2distCoeffs[1,1],k2=C2distCoeffs[1,2],k3=C2distCoeffs[1,5])
    Camera2$Intrinsic$tangential_distortion=list(p1=C2distCoeffs[1,3],p2=C2distCoeffs[1,4])
    return(list(Camera1=Camera1, Camera2=Camera2))
  }
  else {
    stop('Not a valid file type!')
  }

}
