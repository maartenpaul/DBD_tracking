track_msd <- function(x,n=5,framerate=30,pxsize=100){
  x$X <- (x$X*pxsize)/1000
  x$Y <- (x$Y*pxsize)/1000

  coef <- ddply(x,.variables = "track",.fun= function(x) {
    result <- vector(length = n)
    #put first frame of track at zero
    x$frame  <- x$frame-x$frame[1]+1
    if(nrow(x)>n+1){
      sapply(1:n,function(j){


        sum <- unlist(sapply(1:nrow(x),simplify = F,function(k){
          if (is.element(x[k,1]+j,x[,1])){
            which_point <- which(x[,1]==x[k,1]+j)
            sum <- ((x[k,2]-x[which_point,2])^2+(x[k,3]-x[which_point,3])^2)
            return(sum)

        }}))

        result[j] <<- mean(sum)

    })
      return(result)

    }
  })
  return(coef)
}


TRACK_MSD <- function(x,n=5,framerate=30,pxsize=100){
  UseMethod("TRACK_MSD")
}

TRACK_MSD.default <- function(x,n=5,framerate=30,pxsize=100){
  stop("MSD requires data frame")
}

TRACK_MSD.data.frame <-  function(x,n=5,framerate=30,pxsize=100){
  track_msd(x,n,framerate,pxsize)

}

TRACK_MSD.list <-  function(x,n=5,framerate=30,pxsize=100){

  llply(x,function(x){
    out <- TRACK_MSD(x,n,framerate,pxsize)
    out$cellID <- x$cellID
    return(out)
  })

}

