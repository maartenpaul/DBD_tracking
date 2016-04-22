##notes
#interesting statistics
#- MSD points, MSD fit
#- mean speed
# filter for longer tracks


library(plyr)

##Functions
MSD <- function(tracks,n=5,framerate=30,pxsize=100){
  tracks[,2] <- (tracks[,2]*pxsize)/1000
  tracks[,3] <- (tracks[,3]*pxsize)/1000
  
  coef <- ddply(tracks,.variables = "track",.fun= function(x) {
    result <- vector(length = n) 
  #put first frame of track at zero
    x[,1]  <- x[,1]-x[1,1]+1
   if(nrow(x)>n+1){
    for (j in 1:n){
      sum <- 0
      n_sum <- 0
      for(k in 1:(nrow(x))){
        if (is.element(x[k,1]+j,x[,1])){
          which_point <- which(x[,1]==x[k,1]+j)
          sum <- sum+((x[k,2]-x[which_point,2])^2+(x[k,3]-x[which_point,3])^2)
          n_sum <- n_sum+1
        }
      }
    result[j] <- sum/n_sum
   
  }
  return(result)
   }
    
  })
  return(coef)
}

MSD_fit <- function(x,n=5,fitzero=T,framerate=1/33,pxsize=100,offset=0.01){

    #apply function for every track (column V4)
    coef <- adply(x,.margins=1, .expand=FALSE,function(x) {
        mst <- t(x[,-1])
        time <- seq(1,length(mst),1)
        time <- time/framerate/1000 #covert frame to second
        #plot(mst,ylim=c(0,0.5),xlim=c(0,10))

      if(!fitzero){
        mst <- mst-offset
      }
          #fit <- lm(mst[,1] ~ time+0)
          # fit <- tryCatch(nls(formula = mst ~ time*x1+0,lower = list(x1=0,x2=-1),start=list(x1=0,x2=0),algorithm = "port"), error=function(e) NULL)
          fit <- tryCatch(lm(formula = mst ~ time-1), error=function(e) NULL)
          if (!is.null(fit)){
            RSS.p <- sum(residuals(fit)^2)
            TSS <- sum((mst - mean(mst,na.rm = T))^2,na.rm = T)
            rsq <- 1-(RSS.p/TSS)
            
            coef <- c(coef(fit),0.02,rsq)
          } 
         
#           else if (fitzero==FALSE) {
#       #fit <- lm(mst[,1] ~ time)
#       fit <- tryCatch(nls(formula = mst[,1] ~ time*x1+x2,lower = list(x1=0,x2=0),upper = list(x1=100,x2=0.05), start=list(x1=0,x2=0),algorithm = "port"), error=function(e) NULL)
#       
#       if (!is.null(fit)){
#         RSS.p <- sum(residuals(fit)^2)
#         TSS <- sum((mst[,1] - mean(mst[,1]))^2)
#         rsq <- 1-(RSS.p/TSS)
#           coef <- c(coef(fit),rsq)
#       } else {
#         coef <- c(NA,NA,NA)
#       }
#     
#     }
        # MSD=4*D*t
    coef[1] <- coef[1]/4
    #coef[1] <- as.numeric(coef[1]/4)
    return(coef)
     
      
  })
   names(coef) <-  c("track","D","intercept","Rsquared")#,paste("dx_",1:n,sep=""),paste("N_",1:n,sep=""))
  
  return (coef)
}

mean_speed <- function(x,framerate=30,pxsize=100){
  #speed um/ms
  
  #scale to um
  x[,2] <- (x[,2]*pxsize)/1000
  x[,3] <- (x[,3]*pxsize)/1000
  
  #apply function for every track (column V4)
  meansp <- ddply(x,.variables = "track",.fun= function(x) {
  speed <- 0
    for (i in 2:nrow(x)){
        speed <- speed + (x$X[i]-x$X[i-1])^2+(x$Y[i]-x$Y[i-1])^2/((x$frame[i]-x$frame[i-1])*framerate)
      }
    speed <- speed/nrow(x)
  })
  names(meansp) <- c("trackid","speed")
  return(meansp)
}

segment_statistics <- function(segments){
  
  result <- dlply(segments,.variables = c("cellID","track"),function(x){
  seg_angle <- vector()
  for(i in 1:(nrow(x)-2)){
    A <- as.numeric(x[i,2:3])
    B <- as.numeric(x[i+1,2:3])
    C <- as.numeric(x[i+2,2:3])
    AB <- B-A
    CB <- C-B
    dAB <- sqrt((B[1]-A[1])^2+(B[2]-A[2])^2)
    dBC <- sqrt((C[1]-B[1])^2+(C[2]-B[2])^2)
    
    seg_angle <- c(seg_angle,(acos((AB%*%CB)/(dAB*dBC))*180/pi))
    
    
  }
    seg_angle <- c(NA,seg_angle,NA)
    return(seg_angle)
  })
}

chull_stat <- function(x){
  library(pracma)
  x[,2] <- (x[,2]*pxsize)/1000
  x[,3] <- (x[,3]*pxsize)/1000
  out <- ddply(x,.variables = "track",.fun= function(x) {
    x <- cbind(x$X,x$Y)
    y <- chull(x)
    area <- polyarea(x[rev(y),1], x[rev(y),2])
    perimeter <- poly_length(x[rev(y),1], x[rev(y),2])
  D <- princomp(x)
  angle <- atan2(D$loadings[2,1],D$loadings[1,1])
  #return(data.frame("sd"=((sd(x$X)+sd(x$Y))/2)*2.35,"N"=nrow(x),"channel"=1))
  # return(data.frame(,"N"=nrow(x),"channel"=1))
  return(data.frame("N"=nrow(x),"meanX"=mean(x[,1]),"meanY"=mean(x[,2]), "sd"=((sd(x[,1])+sd(x[,2]))/2),"sdpri"=((D$sdev[1]+D$sdev[2])/2),"width"=(max(D$scores[,1])-min(D$scores[,1])),"ratio"=(D$sdev[1]/D$sdev[2]),"angle"=angle,"chull_area"=area,"chull_perimeter"=perimeter))
})
  
return(out)
}

##parameters
#directory <- "D:/Stack/experimental data/Tracking all_data/"
directory <- "D:/Stack/experimental data/160325 RA54-SNAP tracking/200ms tracking"
directory <- "D:/Stack/experimental data/160415 tracking beads_SiR/"
framerate <- 1/200 #ms 1/200
n <- 4 #number of timepoints taken into account for MSD fit
fitzero <- FALSE #should fit go through origin (0,0)
min_length <- 20
pixelsize <- 100 #nm
fitMSD <- F
offset <- 4*0.01^2 #33 nm led to a constant offset (r2(0)) in MSD of 0.0044 mm2



###
condition_list <- list.dirs(directory,full.names = F,recursive = F)
#condition_list <- "SNAP-SiR"

for (i in 1:length(condition_list)){
  segments <- list()
  dir <- file.path(directory,condition_list[i])
  filelist <- list.dirs(dir,full.names = T,recursive = F)
  total <- length(filelist)
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  D <- list()
  
  for(j in 1:total){
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, j)
   # tracks_simple <- 005 SNAP-SiR tracking","tracks.simple.filtered.txt"),sep = "\t",header = F)
    tracks_simple <- read.csv(file.path(filelist[j],"tracks.simple.filtered.txt"),sep = "\t",header = F)
    names(tracks_simple) <- c("frame","X","Y","track","displacement","intensity","sigma","fit_error")
    segments[[basename(filelist[j])]] <- data.frame(tracks_simple,"cellID"=j)
    msd <- MSD(tracks = tracks_simple,n = n,framerate=framerate)
    if(fitMSD){
      coef <-  MSD_fit(msd,n = n,fitzero = fitzero,framerate=framerate,offset)
    #result <- mean_speed(x = tracks_simple,framerate=framerate)
      stats <- chull_stat(x=tracks_simple)
      D[[basename(filelist[j])]] <- data.frame("cellID"=j,msd,coef)#,"speed"=result$speed,stats[,-1])
    }
    D[[basename(filelist[j])]] <- data.frame("cellID"=j,msd)
    write.csv(coef,paste(filelist[j],"D.txt",sep="_"))
  }
  segments <- ldply(segments)
  write.csv(segments,file = file.path(dir,"segments.txt"))
  save(segments,file=file.path(dir,"segments.R"))
  #load(file=file.path(dir,"segments.R"))
  D <- ldply(D)
  D$condition <- condition_list[i]
  save(D,file=file.path(dir,"D.R"))
  close(pb)
}

condlist <- list.dirs(directory,full.names = T,recursive = F)
for (i in 1:length(condlist)){
  load(file.path(condlist[i],"D.R"))
  load(file.path(condlist[[i]],"segments.R"))
  .GlobalEnv[[paste("D_",condition_list[i],sep="")]] <- D
  .GlobalEnv[[paste("segments_",condition_list[i],sep="")]] <- segments
}

allD <- rbind(D_WT_noIR,D_WT_IR,D_KA_noIR,D_KA_IR)
allD <- subset(allD,allD$D>0)
m <- ggplot(allD, aes(x=D,y=..density..,color=condition)) + xlim(c(0.000001,1E3))
m + geom_density() + scale_x_log10() 