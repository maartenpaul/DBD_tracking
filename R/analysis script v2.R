#required packages
library(plyr)
library(lattice)
library(stats)
library(MSDtracking)
library(ggplot2)

# #directory <- "D:/Stack/Genetics/170504 tracking BRCA2/"
# directory <- "D:/Stack/Genetics/170313 Halo dDBD and CTD particle tracking/"
# directory <- "D:/Stack/Genetics/170623 BRCA2 MMC tracking/"
# directory <- "D:/Stack/Genetics/170705 BRCA2 HaloTag Fucci G1_G2/"
# directory <- "D:/Stack/Genetics/170825 BRCA2 PALB2 Halo transfection lipofectamine/"
# directory <- "D:/Stack/Genetics/170602 BRCA2 PALB2 tracking/"
# directory <- "D:/Stack/Genetics/170711 BRCA2 dDBD CTD Halo HU/"
# directory <- "D:/Stack/experimental data/160325 RA54-SNAP tracking/200ms tracking"
# directory <- "D:/Stack/experimental data/160415 tracking beads_SiR/"
# directory <- "D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD/"
#directory <- "D:/Stack/Genetics/171219 Halo BRCA2 alpha tracks/"
directory <- "D:/Stack/Genetics/180110 BRCA2 WTG10_dDdC_F4/"


#input variables
framerate <- 1/32 #1/200 #1/ms
n <- 4 #number of timepoints taken into account for MSD fit
fitzero <- TRUE #should fit go through origin (0,0)

min_length <- 6 #minimum length track
pixelsize <- 100 #nm
fitMSD <- T
offset <- 0.003158194#experimentally determined
max_tracks <- 30 #maximum number of tracks per frame else exclude tracks from dataset, avoids mislinking of tracks

condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- condition_list[c(5)]

segments_all <- list()
msd_fit_all <- list()
track_stats_all <- list()

#load data
for (i in 1:length(condition_list)){
  segments <- list()
  dir <- file.path(directory,condition_list[i])
  filelist <- list.dirs(dir,full.names = T,recursive = F)
  total <- length(filelist)
  # create progress bar
  for(j in 1:total){
  #  Sys.sleep(0.1)
    tracks_simple <- read.csv(file.path(filelist[j],"tracks.simple.filtered.txt"),sep = "\t",header = F)
    names(tracks_simple) <- c("frame","X","Y","track","displacement","intensity","sigma","fit_error")
    segments[[j]] <- data.frame(SEGMENT_STAT(tracks_simple),"cellID"=basename(dirname(file.path(filelist[j],"tracks.simple.filtered.txt"))))
  }
  track_msd <- TRACK_MSD(segments,n = n,framerate=framerate)
  save(track_msd,file=file.path(dir,"track_msd.Rdata"))



  if(fitMSD){
    tracks <-  TRACK_MSD_fit(track_msd,n = n,fitzero = fitzero,framerate=framerate,offset)
    save(tracks,file=file.path(dir,"msd_fit.Rdata"))

    stats <- TRACK_STAT(x=segments)
    save(stats,file=file.path(dir,"track_stats.Rdata"))

  }
  segments_all[[basename(dir)]] <- data.frame(ldply(segments),"condition"=basename(dir))
  msd_fit_all[[basename(dir)]] <- data.frame(ldply(tracks),"condition"=basename(dir))
  track_stats_all[[basename(dir)]] <- data.frame(ldply(stats),"condition"=basename(dir))


}
#save data to the folder
save(segments_all,file=file.path(directory,"segments_all.Rdata"))
save(msd_fit_all,file=file.path(directory,"msd_fit_all.Rdata"))
save(track_stats_all,file=file.path(directory,"track_stats_all.Rdata"))

