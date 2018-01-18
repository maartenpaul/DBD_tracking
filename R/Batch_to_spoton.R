directory <- "D:/Stack/Genetics/171024 PALB2-Halo tracking/"
directory <- "D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD/"

framerate <- 1/30 #1/200 #ms
pixelsize <- 100 #nm

condition_list <- list.dirs(directory,full.names = F,recursive = F)
condition_list <- condition_list[c(4)]
#condition_list <- "SNAP-SiR"
segments_all <- list()
msd_fit_all <- list()
track_stats_all <- list()
#load data
for (i in 1:length(condition_list)){
  segments <- list()
  dir <- file.path(directory,condition_list[i])
  filelist <- list.dirs(dir,full.names = T,recursive = F)
  #filelist <- filelist[-grep("skip",x = filelist)]
  total <- length(filelist)
  # create progress bar
  for(j in 1:total){
    #  Sys.sleep(0.1)
    tracks_simple <- read.csv(file.path(filelist[j],"tracks.simple.filtered.txt"),sep = "\t",header = F)
    spoton <- tracks_simple[c(1,1,4,2,3)]
    spoton[,2] <- spoton[,2]/framerate/1000 #seconds
    spoton[,4:5] <- spoton[,4:5]*pixelsize/1000 #um
    spoton <- cbind(seq(0,nrow(spoton)-1),spoton)
    names(spoton) <- c("","frame","t","trajectory","x","y")
    write.csv(spoton,file = file.path(dirname(filelist[j]),paste(basename(dirname(filelist[j])),j,"tracks.spoton.txt",sep = "_")),row.names=FALSE,quote = FALSE)
  }
}

