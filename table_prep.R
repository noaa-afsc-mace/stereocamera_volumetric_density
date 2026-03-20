rm(list = ls())

#jellyfish
targets<-read.csv("//akc0ss-n086/MACE_Acoustic2/DY2408/Pelagicam/D20240717-T193700/data/DY2308_minicam_pilot_D20240717-T193700_targets.csv",header=TRUE)
deps=unique(targets$DEPLOYMENT_ID)

all_frame_df=data.frame(DEPLOYMENT_ID=character(), FRAME_NUMBER=integer())

for (i in 1:length(deps)){
  ind=which(targets$DEPLOYMENT_ID==deps[i])
  subfr=targets[ind,]
  all_frames=seq(min(subfr$FRAME_NUMBER),max(subfr$FRAME_NUMBER),by=1)
  for (j in 1:length(all_frames)){
    all_frame_df[nrow(all_frame_df) + 1,]=c(deps[i],all_frames[j])
  }

}

counts=vector(mode="numeric", nrow(all_frame_df))
for (i in 1:length(counts)){
  ind=which(targets$FRAME_NUMBER==all_frame_df$FRAME_NUMBER[i] & targets$DEPLOYMENT_ID==all_frame_df$DEPLOYMENT_ID[i])
  counts[i]=length(ind)

}

density=counts/0.2
print(paste("maxN = ",max(counts)))
print(paste("mean density = ",round(mean(density),3)))
print(paste("standard deviation = ",round(sd(density),3)))


#############################################################################

# krill
# targets<-read.csv("C:/Users/kresimir.williams/Work/volumetric desity estimation/stereocamera_volumetric_density/data-raw/krill data.csv",header=TRUE)
# deps=unique(targets$deployment_ID)
#
# all_frame_df=data.frame(deployment_ID=character(), frame_number=integer())
#
# for (i in 1:length(deps)){
#   ind=which(targets$deployment_ID==deps[i])
#   subfr=targets[ind,]
#   all_frames=seq(min(subfr$frame_number),max(subfr$frame_number),by=1)
#   for (j in 1:length(all_frames)){
#     all_frame_df[nrow(all_frame_df) + 1,]=c(deps[i],all_frames[j])
#   }
#
# }
#
# counts=vector(mode="numeric", nrow(all_frame_df))
# for (i in 1:length(counts)){
#   ind=which(targets$frame_number==all_frame_df$frame_number[i] & targets$deployment_ID==all_frame_df$deployment_ID[i])
#   counts[i]=length(ind)
#
# }
#
# density=counts/0.12
# print(paste("maxN = ",max(counts)))
# print(paste("mean density = ",round(mean(density),3)))
# print(paste("standard deviation = ",round(sd(density),3)))


#########################################################################




#pollock
# targets<-read.csv("C:/Users/kresimir.williams/Work/volumetric desity estimation/stereocamera_volumetric_density/data-raw/targets.csv",header=TRUE)
# #targets=targets |> subset(SPECIES_GROUP=="Adult_pollock")
#
#
# # deps=unique(targets$deployment_ID)
# all_frame_df=data.frame(FRAME_NUMBER=integer())
# all_frames=seq(min(targets$FRAME_NUMBER),max(targets$FRAME_NUMBER),by=1)
# cnt=1
# for (j in 1:length(all_frames)){
#   all_frame_df[nrow(all_frame_df) + 1,]=c(all_frames[j])
# }
#
# targets=targets |> subset(SPECIES_GROUP=="Adult_pollock")
# # for (i in 1:length(deps)){
# #   ind=which(targets$deployment_ID==deps[i])
# #   subfr=targets[ind,]
# #   all_frames=seq(min(subfr$frame_number),max(subfr$frame_number),by=1)
# #   for (j in 1:length(all_frames)){
# #     all_frame_df[nrow(all_frame_df) + 1,]=c(deps[i],all_frames[j])
# #   }
# #
# # }
#
# counts=vector(mode="numeric", nrow(all_frame_df))
# for (i in 1:length(counts)){
#   ind=which(targets$FRAME_NUMBER==all_frame_df$FRAME_NUMBER[i])# & targets$deployment_ID==all_frame_df$deployment_ID[i])
#   counts[i]=length(ind)
#
# }
#
# density=counts/5.89
# print(paste("maxN = ",max(counts)))
# print(paste("mean density = ",round(mean(density),3)))
# print(paste("standard deviation = ",round(sd(density),3)))

################################################################################

# yelloweye
# targets<-read.csv("C:/Users/kresimir.williams/Work/volumetric desity estimation/stereocamera_volumetric_density/data-raw/targets3.csv",header=TRUE)
# targets=targets[which(targets$SPECIES_GROUP=="Yelloweye R." & !is.na(targets$RANGE)),]
# rownames(targets)=1:nrow(targets)
#
# for (i in 1:nrow(targets)){
#   depstr=targets$DEPLOYMENT_ID[i]
#   bits=unlist(strsplit(depstr,":"))
#   targets$DEPLOYMENT_ID[i]=bits[[3]]
# }
#
# all_frames=read.csv("C:/Users/kresimir.williams/Work/trigcam/2019_coop_analysis/KW_analysis/2019_trigcam_frames.csv")
# all_frames_yell=data.frame(DEPLOYMENT_ID=character(),FRAME_NUMBER=integer())
#
# deps=unique(targets$DEPLOYMENT_ID)
# for (i in 1:nrow(all_frames)){
#   if (all_frames$deployment_ID[i] %in% deps){
#     all_frames_yell[nrow(all_frames_yell)+1,]=c(all_frames$deployment_ID[i],all_frames$frame_number[i])
#   }
#
# }
#
# counts=vector(mode="numeric", nrow(all_frames_yell))
# for (i in 1:length(counts)){
#   ind=which(targets$FRAME_NUMBER==all_frames_yell$FRAME_NUMBER[i] & targets$DEPLOYMENT_ID==all_frames_yell$DEPLOYMENT_ID[i])
#   counts[i]=length(ind)
#
# }
#
# density=counts/2.69
# print(paste("maxN = ",max(counts)))
# print(paste("mean density = ",round(mean(density),3)))
# print(paste("standard deviation = ",round(sd(density),3)))
