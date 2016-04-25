

track_msd_fit <- function(x,n=5,fitzero=T,framerate=1/33,pxsize=100,offset=0.01){

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

TRACK_MSD_fit <- function(x,y,prec,ch, px,xlim,ylim,file,output,fit,fast){
  UseMethod("TRACK_MSD_fit")
}

TRACK_MSD_fit.default <- function(x,n=5,fitzero=T,framerate=1/33,pxsize=100,offset=0.01){
  stop("MSD_fit requires data frame")
}

TRACK_MSD_fit.data.frame <-  function(x,n=5,fitzero=T,framerate=1/33,pxsize=100,offset=0.01){
  track_msd_fit(x,n,fitzero,framerate,pxsize,offset)

}

TRACK_MSD_fit.list <-  function(x,n=5,fitzero=T,framerate=1/33,pxsize=100,offset=0.01){
  llply(x,function(x){
    TRACK_MSD_fit(x,n,fitzero,framerate,pxsize,offset)
  })
}

