track_msd <- function(tracks,n=5,framerate=30,pxsize=100){
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


TRACK_MSD <- function(tracks,n=5,framerate=30,pxsize=100){
  UseMethod("TRACK_MSD")
}

TRACK_MSD.default <- function(tracks,n=5,framerate=30,pxsize=100){
  stop("MSD requires data frame")
}

TRACK_MSD.data.frame <-  function(tracks,n=5,framerate=30,pxsize=100){
  track_msd(tracks,n,framerate,pxsize)

}

TRACK_MSD.list <-  function(tracks,n=5,framerate=30,pxsize=100){
  llply(x,function(x){
    TRACK_MSD(tracks,n,framerate,pxsize)
  })
}

