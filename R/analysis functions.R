msd_analyze_data <- function(directory,condition_list,framerate,n,fitzero,min_length,pixelsize,fitMSD,offset,max_tracks){
  segments_all <- list()
  msd_fit_all <- list()
  track_stats_all <- list()

  for (i in 1:length(condition_list)){
    segments <- list()
    dir <- file.path(directory,condition_list[i])
    filelist <- list.dirs(dir,full.names = T,recursive = F)
    total <- length(filelist)
    # create progress bar
    for(j in 1:total){
      Sys.sleep(0.1)
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
}
theme_Publication <- function(base_size=10, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="histdatack"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_histdatank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

msd_histogram <- function(msd_fit_all,directory,name="",threshold=0.05,order=NULL,ymax=0.12){


  tracks <- ldply(msd_fit_all)
  out <- ddply(tracks,.variables = c("condition","cellID")  ,function(x){
    out <- hist(log10(x$D),breaks = seq(-5,2,0.1),plot = F)
    return(out$counts/sum(out$counts))
  })


  means <- ddply(out,.variables = "condition",function(x){
    colMeans(x[,-c(1,2)])
  })
  sds <- ddply(out,.variables = "condition",function(x){
    apply(x[,-c(1,2)], 2, sd)#/sqrt(length(table(x$cellID))) #if se instead of stdev
  })

  means <- melt(means)
  sds <- melt(sds)
  mids <- seq(-4.9,2,0.1)
  mids <- 10^mids

  histdata <- data.frame(".id"=means$condition,"x"=rep(mids,each=length(msd_fit_all)),'mean'=means$value,'se'=sds$value)
  if(!is.null(order)){
    histdata$.id <- factor(histdata$.id,levels(histdata$.id)[order])
  }

  #histdata$se <- 0
  histdata <- na.omit(histdata)

  data <- data.frame(histdata[,1:2])
  names(data) <- c("x","y")

  #levels(histdata$.id) <- c(levels(histdata$.id)[2],levels(histdata$.id)[1])

  q1 <- ggplot(histdata, aes(x=x, y=mean,colour=.id,fill=.id)) +
    geom_smooth(stat="identity")+
    # geom_smooth(stat="smooth",span = 0.2)+
    geom_ribbon(aes(ymin=mean-se, ymax=mean+se,colour=NULL), alpha=0.4)+ylim(0,ymax)+
    theme(text = element_text(size=20),legend.position = "none")+facet_wrap(~.id) +scale_x_continuous(trans="log10")+ scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+theme(legend.position = "none")+ geom_vline(aes(xintercept=threshold),color = "red",linetype="dashed")
  print(q1+xlab(expression(D[eff]~mu~m^{2}/s)))

  ggsave(filename = file.path(directory,paste0(name,"_D_histogram.pdf")),plot = q1)

  results <- llply(msd_fit_all,function(x) {
    daply(x,.variables = "cellID",function(x){
      out <- table(x$D>threshold)/length(x$D)
      return(out)
    })
  }
  )

  results2 <- ldply(results)
  names(results2) <- c(".id","immobile","mobile")
  write.table(results2,file=file.path(directory,paste0(name,"_imm_fractions.txt")),,row.names = F,col.names = T)

  if(!is.null(order)){
    results2$.id <- factor(results2$.id,levels(histdata$.id)[order])
  }
  results2 <- na.omit(results2)
  out <- ldply(results, function(x){
    means <- mean(x[,1])
    sems <- sd(x[,1])/sqrt(length(x[,1]))
    return(data.frame("mean"=means,"se"=sems))
  }
  )
  #out$.id <- factor(out$.id,out$.id[c(1,2)])


  #out$.id <- factor(out$.id,levels(histdata$.id)[c(1,2,3,4)])
  out <- na.omit(out)
  q2 <- ggplot(out, aes(x=.id, y=mean,fill=.id)) +
    geom_bar(position=position_dodge(width = 1.1), stat="identity")+  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2)+xlab("")+ylab("immobile fraction")+
    scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=16)+theme(legend.position = "none")+ylim(0,0.5)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(q2)

  ggsave(filename = file.path(directory,paste0(name,"_immobile_bar_graph.pdf")),plot = q2)

  q3 <- ggplot(results2,aes(x=factor(.id),y=immobile,fill=.id))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_dotplot(binaxis = "y", stackdir = "center")+ scale_colour_Publication()+scale_fill_Publication()+theme(legend.position = "none")+
    ylim(0,0.5)
  ggsave(filename = file.path(directory,paste0(name,"_immobile_dot_plot.pdf")),plot = q3)


  print(q3)
  #qa <- grid.arrange(q1,q2,ncol=2)
}

sos_to_spoton_csv <- function(condition_list,directory,framerate,pixelsize){
  condition_list <- list.dirs(directory,full.names = F,recursive = F)
  #condition_list <- condition_list[c(4)]
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
}
