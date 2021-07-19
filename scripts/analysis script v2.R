#required packages
library(plyr)
library(lattice)
library(stats)
library(MSDtracking)
library(ggplot2)
library(ggpol)

#input variables
framerate <- 1/32 #1/200 #1/ms
n <- 4 #number of timepoints taken into account for MSD fit
fitzero <- TRUE #should fit go through origin (0,0)

min_length <- 6 #minimum length track to be processed
pixelsize <- 100 #nm
fitMSD <- T
offset <- 4*(0.01)^2#experimentally determined
max_tracks <- 500 #maximum number of tracks per frame else exclude tracks from dataset, avoids mislinking of tracks

directory <- ""
condition_list <- list.dirs(directory,full.names = F,recursive = F)

msd_analyze_data(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks)

load(file.path(directory,"msd_fit_all.Rdata"))
load(file.path(directory,"segments_all.Rdata"))

#plot the results from basic MSD analysis
msd_histogram(msd_fit_all,directory,threshold = 0.05)
