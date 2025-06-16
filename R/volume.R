#' Get volume from a Matlab input list. Everything hard-coded still FIX ME
#' @param matcal A list of calibration values from the \code{readMat} function.
#' @return A list containing the function vol_func and a vector of coefficients p_vol
#' @references Bouget (2008) Camera calibration toolbox for Matlab. Available from http://vision.caltech.edu/bouguetj/calib_doc/index.html (accessed September 2008)].
#' @details to do
#' @export
get_vol_func <- function(matcal){
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
  p_vol= as.numeric(coef(lm(vol ~cam_range +I(cam_range^2))))
  disp_range=seq(0,7,length.out = 100)
  y=p_vol[1]+p_vol[2]*disp_range+p_vol[3]*disp_range^2

  # create the function object for itegration
  vol_func<- function(x){p_vol[1]+p_vol[2]*x+p_vol[3]*x^2}
  return(list(vol_func=vol_func, p_vol=p_vol, vol=vol))
}
