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

min_length <- 6 #minimum length track
pixelsize <- 100 #nm
fitMSD <- T
offset <- 4*(0.01)^2#experimentally determined
max_tracks <- 500 #maximum number of tracks per frame else exclude tracks from dataset, avoids mislinking of tracks

#directory <- "D:/Stack/Genetics/180823 BRCA2_Halo dCTD dDBD HU/"
directory <- "F:/181012 BRCA2-Halo PCNA-iRFP IR tracking/"
#directory <- "O:/Maarten/Genetics/180117 BRCA2 dDdC IR tracking"
condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- condition_list[c(1,2)]


msd_analyze_data(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks)

load(file.path(directory,"msd_fit_all.Rdata"))
#load(file.path(directory,"track_stats_all.Rdata"))
load(file.path(directory,"segments_all.Rdata"))

msd_histogram(msd_fit_all,directory,threshold = 0.05)

#cell cycle
msds <- msd_fit_all$`WT G10 cell cycle`

msd_cell_cycle <- list()
msd_cell_cycle$G1 <- subset(msds,cellID%in%c("002 30ms","005 30ms","007 30ms","010 30ms"))
msd_cell_cycle$G1$condition <- "G1"
msd_cell_cycle$S <- subset(msds,cellID%in%c("003 30ms","004 30ms","006 30ms","008 30ms"))
msd_cell_cycle$S$condition <- "S"
msd_cell_cycle$G2 <- subset(msds,cellID%in%c("001 30ms","009 30ms"))
msd_cell_cycle$G2$condition <- "G2"

msd_histogram(msd_cell_cycle,directory,threshold = 0.05)
