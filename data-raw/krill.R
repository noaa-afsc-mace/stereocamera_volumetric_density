## code to prepare 'krill' dataset goes here

if(!requireNamespace("yaml", quietly = TRUE))
  stop("yaml package needed to read in calibration for the krill data set")

# reading in matlab file with stereo calibration parameters
cal <- yaml.load(read_yaml("calibrations/example_calibration_frommatlab2.yml"))

######################################################################
#  TARGET DETECTION FUNCTION ############

# read in krill data
targets<-read.csv("data-raw/krill data.csv",header=TRUE)
targets=targets |> subset(class=="Krill")
targets$Range=targets$Range/10# change units from cm to dm ( this gives volume in L)
krill <- list(cal=cal, targets=targets)

usethis::use_data(krill, overwrite = TRUE)
