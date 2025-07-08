
rm(list = ls())
library(yaml)
readCalibration = function(filename,method='matlab_caltech'){
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
    Camera1$Intrinsic$focal_point=list(h=Calraw$fc.left[1],v=Calraw$fc.left[2])
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
    Camera2$Intrinsic$focal_point=list(h=Calraw$fc.right[1],v=Calraw$fc.right[2])
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
    Camera1$Intrinsic$focal_point=list(h=fl_x,v=fl_y)
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
    Camera2$Intrinsic$focal_point=list(h=fl_x,v=fl_y)
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
    Camera1$Intrinsic$focal_point=list(h=C1Matrix[1,1],v=C1Matrix[2,2])
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
    Camera2$Intrinsic$focal_point=list(h=C2Matrix[1,1],v=C2Matrix[2,2])
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


Cal=readCalibration('accessory_functions/example_calibrations/matlab_caltech_example.mat','matlab_caltech')
z <- as.yaml(list(Camera1=Cal$Camera1,Camera2=Cal$Camera2))
write_yaml(z,'calibrations/example_calibration_frommatlab.yml')

Cal=readCalibration(list('accessory_functions/example_calibrations/seagis_cal_example_Left.TXT','accessory_functions/example_calibrations/seagis_cal_example_Right.TXT'),'seagis')
z <- as.yaml(list(Camera1=Cal$Camera1,Camera2=Cal$Camera2))
write_yaml(z,'calibrations/example_calibration_fromseagis.yml')

# Cal=readCalibration('accessory_functions/example_calibrations/opencv_example.npz','opencv')
# z <- as.yaml(list(Camera1=Cal$Camera1,Camera2=Cal$Camera2))
# write_yaml(z,'calibrations/example_calibration_fromopencv.yml')
