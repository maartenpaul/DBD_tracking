#required packages
library(plyr)
library(lattice)
library(stats)
library(MSDtracking)
library(ggplot2)

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
condition_list <- condition_list[c(5)]

directory <- "D:/Stack/Genetics/180117 BRCA2 dDdC IR tracking/"

msd_analyze_data(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks)

msd_histogram(msd_fit_all,directory)
