library(tidyverse)
library(plyr)
library(doParallel)

library(ggpol)
theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
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
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#ffffff"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.4, "cm"),
            legend.margin = unit(0.2, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#ffffff",fill="#ffffff"),
            strip.text = element_text(face="bold")
    ))

}


scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#c00000","#fdae61","#1f497d","#6599d9","#542788","#de77ae","#271d68","#6dc5aa")), ...)

}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#c00000","#fdae61","#1f497d","#6599d9","#542788","#de77ae","#217d68","#6dc5aa")), ...)

}

#run ML MSS
datasets <- c("/media/OIC-station2/Genetics/181005 BRCA2-Halo PCNA-iRFP IR tracking/","/media/OIC-station2/Genetics/181010 BRCA2-Halo PCNA-iRFP IR tracking/","/media/OIC-station2/Genetics/181012 BRCA2-Halo PCNA-iRFP IR tracking/")
segments <- list()
for(x in datasets){
  load(file.path(x,"segments_all.Rdata"))
  segments_all <- ldply(segments_all,.id = NA)
  segments_all$.id <- NULL
  segments[[basename(x)]] <- segments_all
}

segments_all <- ldply(segments)
rm(segments)
segments_all$cellID <- paste0(segments_all$cellID,segments_all$.id)
#calculate MSD and MSS


source('python/ML_py.R')

segments_all <- ddply(segments_all, .variables = c("condition","cellID"), function(x){
    ML_segment_tracks(x)
})
save(segments_all,file = "/home/maarten/Documents/DBD paper analysis/segments_all.Rdata")

ptm <- proc.time()
#initialize cluster
nodes <- detectCores()
cl <- makeCluster(nodes-12)
registerDoParallel(cl)

#make tracklet column
  segments_all <- ddply(segments_all,.variables = c("condition","cellID"),.parallel = T,function(x){
    ddply(x,.variables = c("track"), .parallel = F, function(x){

      segment <- 1
      tracklets <- vector(length=nrow(x))
      tracklets[1] <- segment

      for (i in 2:nrow(x)){
        if(x$state[i]==x$state[i-1]){
          tracklets[i] <-segment
        } else {
          segment <- segment+1
          tracklets[i] <-segment
        }
      }
      x$tracklet <- paste0(x$track,".",tracklets)

      return(x)

    })})


stopCluster(cl)
proc.time() - ptm

numPmsd <- 4
numPmss <- 4
minLen <- 6
p <- seq(from=0.5,to=6,length.out=12)
py$pixSize <- 0.100
py$t <- 0.032
source_python('python/getMSDandMSS_R.py')

MSD_MSS <- function(x){
  if(nrow(x)>10){
    out <- getMSDandMSS_R(x$X,x$Y)
    return(tibble("D_ML"=out[[1]],"D_Smmss"=out[[2]]))
  } else {
    return(tibble("D_ML"=-1.0,"D_Smmss"=-1.0))
    return(NA)
  }
}

MSD_only <- function(x){
  if(nrow(x)>10){
    out <- getMSDandMSS_R(x$X,x$Y)
    return(out[[1]])
  } else {
    return(NA)
  }
}

library(tidyverse)
segs_nest <- segments_all %>%
  select(cellID,condition,tracklet,X,Y) %>%
  group_by(cellID,condition,tracklet) %>%
  group_modify(~MSD_MSS(.x)) %>%
  inner_join(y=segments_all,by=c("cellID","condition","tracklet"))

save(segs_nest,file = "/home/maarten/Documents/DBD paper analysis/segments_tbl.Rdata")
write_csv(segs_nest,path = "/home/maarten/Documents/DBD paper analysis/segments_tbl.txt")

segs_nest <- segs_nest %>%
  filter(condition!="WT G10 cell cycle")

segs_nest$condition <- factor(segs_nest$condition,levels=c("WT G10 -IR","WT G10 +IR","dDBD E4 -IR","dDBD E4 +IR","dCTD A2 -IR","dCTD A2 +IR","dDBDdCTD F4 -IR","dDBDdCTD F4 +IR"))

#showing the diffusion rate for fast fractions
p <- segs_nest %>%
  filter(D_ML>0,state==0)%>%
  filter(D_ML>0,grepl("-IR",condition))%>%
    dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
    group_by(condition,cellID)%>%
    dplyr::summarise(mean=mean(D_ML)) %>%
  group_by(condition)%>%
  ggplot(aes(x=condition,y=mean, fill=condition))+geom_boxplot()+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=18)+ theme(legend.position = "none")+
  xlab("")+ylab(expression(D[app]~mu~m^{2}/s))+ylim(0,3)
p
ggsave(p,filename = "plots/diffusionrate_fast.pdf",width = 8,height = 8,units = "in" )


#only -IR
segs_nest %>%
  filter(D_ML>0,grepl("-IR",condition))%>%
  dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
  group_by(condition,cellID,state)%>%
  dplyr::summarise(mean=mean(D_ML)) %>%
  group_by(condition)%>%
  ggplot(aes(x=condition,y=mean,fill=condition))+geom_boxplot(notch = F)+facet_wrap(.~state,scales = "free_y")+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+ theme(legend.position = "none") + ylab(expression(D[app]~mu~m^{2}/s)) + xlab("")

fractions <- segs_nest %>%
  filter(D_ML>0)%>%
  dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
  group_by(condition,cellID,state)%>%
  dplyr::summarise(number=n())%>%
  group_by(condition,cellID) %>%
  dplyr::mutate(fraction=number/sum(number))

plt <- fractions %>%
  filter(state==2) %>%
ggplot(aes(y=fraction,fill=condition))+geom_boxjitter(errorbar.draw = TRUE,jitter.height = 0, jitter.width = 0.015)+ylab("immobile fraction")+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)
plt
ggsave(plt,filename = "plots/immobile_fractions.pdf",width = 8,height = 8,units = "in" )




fractions %>%
  filter(state==2) %>%
  mutate(expID=str_split(fractions$cellID,pattern = "30m",simplify = T,)[3])%>%
  ggplot(aes(y=fraction,x=condition,fill=condition))+geom_dotplot(binaxis='y', stackdir='center')+ylab("immobile fraction")+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+ stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
                                                                                                    geom="errorbar", color="black", width=0.5) +
  stat_summary(fun.y=mean, geom="point", color="black")


x <- fractions %>%
  filter(state==2)

t.test(subset(x,condition=="dD -IR")$fraction,subset(x,condition=="WT G10 +IR")$fraction)

  group_by(condition)%>%
  dplyr::summarise(median=median(mean),sd(mean))%>%
  ggplot(aes(x=median,y=condition))+geom_col()


k <- segs_nest %>%
  filter(D_ML>0)%>%
  dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
  mutate(state_str=as.character(state))

x <- group_by(k,condition) %>%
  group_by(condition,state)%>%
  dplyr::summarise(number=n())%>%
  group_by(condition) %>%
  dplyr::mutate(fraction=number/sum(number))

y <- group_by(k,condition) %>%
  group_by(condition,cellID,state)%>%
  dplyr::summarise(number=n())%>%
  group_by(condition,cellID) %>%
  dplyr::mutate(fraction=number/sum(number)) %>%
  group_by(condition,state) %>%
  dplyr::summarise(mean_fraction=round(mean(fraction),digits = 2),sd_fraction=round(sd(fraction),digits = 2))



p <- segs_nest %>%
  filter(D_ML>0)%>%
  dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
    mutate(state_str=as.character(state))%>%
  ggplot(aes(x=D_ML,fill=state_str,y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+geom_histogram(position="identity",alpha=0.5)+scale_x_log10(limits=c(0.0001,10))+facet_wrap(.~condition,ncol=2)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+ theme(legend.position = "none")+ylab("relative frequency")+ xlab(expression(D[app]~mu~m^{2}/s)) +
  geom_text(data=subset(y,state==0 ), aes(x=1., y=0.12, label=mean_fraction), colour="#c00000", inherit.aes=FALSE, parse=FALSE)+
  geom_text(data=subset(y,state==0 ), aes(x=1., y=0.10, label=paste0("+/-",sd_fraction)), colour="#c00000", inherit.aes=FALSE, parse=FALSE)+

  geom_text(data=subset(y,state==1 ), aes(x=0.08, y=0.06, label=mean_fraction), colour="#fdae61", inherit.aes=FALSE, parse=FALSE)+
  geom_text(data=subset(y,state==1), aes(x=.08, y=0.04, label=paste0("+/-",sd_fraction)), colour="#fdae61", inherit.aes=FALSE, parse=FALSE)+
  geom_text(data=subset(y,state==2 ), aes(x=0.003, y=0.09, label=mean_fraction), colour="#1f497d", inherit.aes=FALSE, parse=FALSE)+
geom_text(data=subset(y,state==2 ), aes(x=.003, y=0.07, label=paste0("+/-",sd_fraction)), colour="#1f497d", inherit.aes=FALSE, parse=FALSE)

p
ggsave(p,filename = "plots/diffusionhistograms.pdf",width = 8,height = 8,units = "in" )
'
p <- segs_nest %>%
  filter(D_Smmss>0)%>%
  dplyr::distinct(condition,cellID,tracklet,.keep_all=T)%>%
  mutate(state_str=as.character(state))%>%
  ggplot(aes(x=D_Smmss,fill=state_str,y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]))+geom_histogram(position="identity",alpha=0.5)+scale_x_log10(limits=c(0.0001,10))+facet_wrap(.~condition,ncol=2)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+ theme(legend.position = "none")+ylab("relative frequency")+ xlab(expression(D[app]~mu~m^{2}/s))
p


#apparent dwell times
segs_nest %>%
  filter(D_ML>0) %>%
  group_by(condition,cellID,tracklet,state)%>%
  dplyr::summarize(length=n())%>%
  mutate(state_str=as.character(state))%>%
  ggplot(aes(x=length,fill=state_str))+geom_density(position="identity",alpha=0.5)+facet_wrap(.~condition,ncol=2)+
  scale_colour_Publication()+scale_fill_Publication()+theme_Publication(base_size=12)+ theme(legend.position = "none")+ylab("relative frequency")+ xlab(expression(D[app]~mu~m^{2}/s))

x<- segs_nest %>%
  filter(D_ML>0) %>%
  group_by(condition,cellID,tracklet,state)%>%
  dplyr::summarize(length=n())%>%
  group_by(condition,state)%>%
  dplyr::summarize(dwelltime=mean(length))
