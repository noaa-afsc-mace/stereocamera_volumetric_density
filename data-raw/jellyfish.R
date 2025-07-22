## code to prepare 'jellyfish' dataset goes here

if(!requireNamespace("yaml", quietly = TRUE))
  stop("yaml package needed to read in calibration for the jellyfish data set")
library(yaml)
# reading in matlab file with stereo calibration parameters
cal <- yaml.load(read_yaml("data-raw/example_calibration_frommatlab.yml"))

######################################################################
#  TARGET DETECTION FUNCTION ############

# read in fish data
targets<-read.csv("data-raw/targets2.csv",header=TRUE)
targets=targets |> subset(SPECIES_GROUP=="Sarsia sp")
targets$RANGE=targets$RANGE/100# change units from cm to m
jellyfish <- list(cal=cal, targets=targets)

usethis::use_data(jellyfish, overwrite = TRUE)
