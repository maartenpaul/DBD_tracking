track_mean_speed <- function(x,framerate=30,pxsize=100){
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

TRACK_MEAN_SPEED <- function(tracks,framerate=30,pxsize=100){
  UseMethod("TRACK_MEAN_SPEED")
}

TRACK_MEAN_SPEED.default <- function(tracks,framerate=30,pxsize=100){
  stop("TRACK_MEAN_SPEED requires data frame")
}

TRACK_MEAN_SPEED.data.frame <-  function(tracks,framerate=30,pxsize=100){
  track_mean_speed(tracks,framerate,pxsize)

}

TRACK_MEAN_SPEED.list <-  function(tracks,framerate=30,pxsize=100){
  llply(x,function(x){
    TRACK_MEAN_SPEED(tracks,framerate,pxsize)
  })
}


track_stat <- function(tracks,framerate=30,pxsize=100){
  library(pracma)
  tracks$X <- (tracks$X*pxsize)/1000
  tracks$Y <- (tracks$Y*pxsize)/1000
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

TRACK_STAT <- function(tracks,framerate=30,pxsize=100){
  UseMethod("TRACK_STAT")
}

TRACK_STAT.default <- function(tracks,framerate=30,pxsize=100){
  stop("TRACK_STAT requires data frame")
}

TRACK_STAT.data.frame <-  function(tracks,framerate=30,pxsize=100){
  track_mean_speed(tracks,framerate,pxsize)

}

TRACK_STAT.list <-  function(tracks,framerate=30,pxsize=100){
  llply(x,function(x){
    TRACK_STAT(tracks,framerate,pxsize)
  })
}


segment_stat <- function(tracks){

  get_angle <- function(x){
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
    seg_angle <- c(-1,seg_angle,-1)
    return(seg_angle)
  }

  result <- ddply(tracks,.variables = "track",function(x){ x$angle <- get_angle(x)
  return(x)
  }
  )
}

SEGMENT_STAT.default <- function(tracks){
  stop("SEGMENT_STAT requires data frame")
}

SEGMENT_STAT.data.frame <-  function(tracks){
  segment_stat(tracks)

}

SEGMENT_STAT.list <-  function(tracks){
  llply(x,function(x){
    SEGMENT_STAT(tracks)
  })
}
