
rm(list = ls())
library(yaml)

# Cal=readCalibration('accessory_functions/example_calibrations/matlab_caltech_example2.mat','matlab_caltech')
# z <- as.yaml(list(Camera1=Cal$Camera1,Camera2=Cal$Camera2))
# write_yaml(z,'calibrations/example_calibration_frommatlab2.yml')

# Cal=readCalibration(list('accessory_functions/example_calibrations/seagis_cal_example_Left.TXT','accessory_functions/example_calibrations/seagis_cal_example_Right.TXT'),'seagis')
# z <- as.yaml(list(Camera1=Cal$Camera1,Camera2=Cal$Camera2))
# write_yaml(z,'calibrations/example_calibration_fromseagis.yml')

# Cal=readCalibration('accessory_functions/example_calibrations/opencv_example.npz','opencv')
# z <- as.yaml(list(Camera1=Cal$Camera1,Camera2=Cal$Camera2))
# write_yaml(z,'calibrations/example_calibration_fromopencv.yml')

Cal=read_calibration('C:/Users/kresimir.williams/Work/volumetric desity estimation/matlab/Pelagicam/2023/trigcam_2019_unit_3.mat','matlab_caltech')
z <- as.yaml(list(Camera1=Cal$Camera1,Camera2=Cal$Camera2))
write_yaml(z,'calibrations/example_calibration_frommatlab3.yml')
