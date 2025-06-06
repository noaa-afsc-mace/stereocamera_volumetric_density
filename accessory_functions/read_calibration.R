
rm(list = ls())

readCalibration = function(filename,method='matlab_caltech'){
  if (method=='matlab_caltech'){
    library(R.matlab)
    Cal = readMat(filename)
    return(Cal)
  }

  else if (method=='seagis'){
    
    
  }
}

Cal1=readCalibration('accessory_functions/example_calibrations/matlab_caltech_example.mat','matlab_caltech')

