## code to prepare 'pollock' dataset goes here

library(yaml)

# reading in matlab file with stereo calibration parameters
cal <- yaml.load(read_yaml("calibrations/example_calibration_frommatlab.yml"))

######################################################################
#  TARGET DETECTION FUNCTION ############

# read in fish data
targets<-read.csv("data-raw/targets.csv",header=TRUE)
targets=targets |> subset(SPECIES_GROUP=="Adult_pollock")
targets$RANGE=targets$RANGE/100# change units from cm to m
pollock <- list(cal=cal, targets=targets)

usethis::use_data(pollock, overwrite = TRUE)
