all <- list()

directory <- "D:/Stack/Genetics/171208 N-BRCA2 Halo H8 MMC/"
load(file.path(directory,"msd_fit_all.Rdata"))
msd_fit_all$`-MMC`$cellID <- paste0('171208_',msd_fit_all$`-MMC`$cellID)
msd_fit_all$`+MMC`$cellID <- paste0('171208_',msd_fit_all$`+MMC`$cellID)

all$`-MMC` <- msd_fit_all$`-MMC`
all$`+MMC` <- msd_fit_all$`+MMC`

directory <- "D:/Stack/Genetics/171201 N-BRCA2 Halo H8 MMC/"
load(file.path(directory,"msd_fit_all.Rdata"))
msd_fit_all$`-MMC`$cellID <- paste0('171201_',msd_fit_all$`-MMC`$cellID)
msd_fit_all$`+MMC`$cellID <- paste0('171201_',msd_fit_all$`+MMC`$cellID)

all$`-MMC` <- rbind(all$`-MMC`,msd_fit_all$`-MMC`)
all$`+MMC` <- rbind(all$`+MMC`,msd_fit_all$`+MMC`)

directory <- "D:/Stack/Genetics/171117 N-BRCA2 MMC/"
load(file.path(directory,"msd_fit_all.Rdata"))
msd_fit_all$`H8 -MCC`$cellID <- paste0('171117_',msd_fit_all$`H8 -MCC`$cellID)
msd_fit_all$`H8 +MCC`$cellID <- paste0('171117_',msd_fit_all$`H8 +MCC`$cellID)

all$`-MMC` <- rbind(all$`-MMC`,msd_fit_all$`H8 -MCC`)
all$`+MMC` <- rbind(all$`+MMC`,msd_fit_all$`H8 +MCC`)

all$`-MMC`$condition <- "-MMC"
all$`+MMC`$condition <- "+MMC"
msd_fit_all <-all

file.path("D:/Stack/Genetics/Tracking data merged/")
