## code to prepare 'pollock' dataset goes here

if(!requireNamespace("yaml", quietly = TRUE))
  stop("yaml package needed to read in calibration for the pollock data set")

# reading in matlab file with stereo calibration parameters
cal <- yaml.load(read_yaml("data-raw/example_calibration_frommatlab.yml"))

######################################################################
#  TARGET DETECTION FUNCTION ############

# read in fish data
targets<-read.csv("data-raw/targets.csv",header=TRUE)
targets=targets |> subset(SPECIES_GROUP=="Adult_pollock")
targets$RANGE=targets$RANGE/100# change units from cm to m
pollock <- list(cal=cal, targets=targets)

usethis::use_data(pollock, overwrite = TRUE)
