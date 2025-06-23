## code to prepare 'pollock' dataset goes here

if(!requireNamespace("R.matlab", quietly = TRUE))
  stop("R.matlab package needed to preapre the pollock data set")

# reading in matlab file with stereo calibration parameters
matcal <- R.matlab::readMat("data-raw/minicam_09262023_cal_sebastes.mat")

######################################################################
#  TARGET DETECTION FUNCTION ############

# read in fish data
targets<-read.csv("data-raw/targets.csv",header=TRUE)
targets=targets |> subset(SPECIES_GROUP=="Adult_pollock")
targets$RANGE=targets$RANGE/100# change units from cm to m
pollock <- list(matcal=matcal, targets=targets)

usethis::use_data(pollock, overwrite = TRUE)
