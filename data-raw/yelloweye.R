## code to prepare 'yelloweye' dataset goes here

if(!requireNamespace("yaml", quietly = TRUE))
  stop("yaml package needed to read in calibration for the pollock data set")

# reading in matlab file with stereo calibration parameters
cal <- yaml.load(read_yaml("calibrations/example_calibration_frommatlab3.yml"))

######################################################################
#  TARGET DETECTION FUNCTION ############

# read in fish data
targets<-read.csv("data-raw/targets3.csv",header=TRUE)
targets=targets[which(targets$SPECIES_GROUP=="Yelloweye R." & !is.na(targets$RANGE)),]
rownames(targets)=NULL
targets$RANGE=targets$RANGE/100# change units from cm to m
yelloweye <- list(cal=cal, targets=targets)

usethis::use_data(yelloweye, overwrite = TRUE)
