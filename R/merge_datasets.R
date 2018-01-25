all <- list()

merge_dataset <- function(dirs,conds,con_name){
  if(length(dirs)!=length(conds)){
    stop("dirs not equal to conds")
  }
  for(i in 1:length(dirs)){
    load(file.path(dirs[i],"msd_fit_all.Rdata"))
    date <- basename(dirs[i])
    date <- strsplit(date," ")
    date <- date[[1]][1]
    msd_fit_all[[conds[i]]]$cellID <- paste0(date,"_",msd_fit_all[[conds[i]]]$cellID)
    if (i==1) {
      all <- msd_fit_all[[conds[i]]]
    } else {
      all <- rbind(all,msd_fit_all[[conds[i]]])
    }
  }
  all$condition <- con_name
  return(all)
}
merge <- merge_dataset(c("D:/Stack/Genetics/180117 BRCA2 dDdC IR tracking",
                "D:/Stack/Genetics/171208 N-BRCA2 Halo H8 MMC"),c("WT G10 +IR","+MMC"),con_name="WT")


###HU
all_data <- list()
all_data[["WTnoHU"]] <- merge_dataset(c("D:/Stack/Genetics/180110 BRCA2 WTG10_dDdC_F4",
                         "Y:/Maarten/Genetics/170711 BRCA2 dDBD CTD Halo HU"),c("WT G10 noHU","WT E10 noHU"),con_name="WT -HU")

all_data[["WTHU"]] <- merge_dataset(c("D:/Stack/Genetics/180110 BRCA2 WTG10_dDdC_F4",
                      "Y:/Maarten/Genetics/170711 BRCA2 dDBD CTD Halo HU"),c("WT G10 +HU 2h","WT E10 HU"),con_name="WT +HU")

all_data[["dDBDdCTDnoHU"]] <- merge_dataset(c("D:/Stack/Genetics/180110 BRCA2 WTG10_dDdC_F4",
                          "Y:/Maarten/Genetics/170711 BRCA2 dDBD CTD Halo HU"),c("dDBDdCTD F4 noHU","dDBDdCTD G9 noHU"),con_name="dDBDdCTD -HU")

all_data[["dDBDdCTDHU"]] <- merge_dataset(c("D:/Stack/Genetics/180110 BRCA2 WTG10_dDdC_F4",
                        "Y:/Maarten/Genetics/170711 BRCA2 dDBD CTD Halo HU"),c("dDBDdCTD F4 +HU 2h","dDBDdCTD G9 HU"),con_name="dDBDdCTD +HU")

all_data[["dDBDdCTDnoHU"]] <- merge_dataset(c("D:/Stack/Genetics/180110 BRCA2 WTG10_dDdC_F4",
                                              "Y:/Maarten/Genetics/170711 BRCA2 dDBD CTD Halo HU"),c("dDBDdCTD F4 noHU","dDBDdCTD G9 noHU"),con_name="dDBDdCTD -HU")

all_data[["dDBDdCTDHU"]] <- merge_dataset(c("D:/Stack/Genetics/180110 BRCA2 WTG10_dDdC_F4",
                                            "Y:/Maarten/Genetics/170711 BRCA2 dDBD CTD Halo HU"),c("dDBDdCTD F4 +HU 2h","dDBDdCTD G9 HU"),con_name="dDBDdCTD +HU")



msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="HU",threshold=0.05)

###MMC
all_data <- list()
datasets <- c("Y:/Maarten/Genetics/170517 BRCA2 tracking MMC","Y:/Maarten/Genetics/170623 BRCA2 MMC tracking","D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD")
all_data[["WTnoMMC"]] <- merge_dataset(datasets,c("WT E10 5nM noMMC","WT E10 noMMC2h","WT E10 noMMC"),con_name="WT -MMC")

all_data[["WTMMC"]] <- merge_dataset(datasets,c("WT E10 5nM 1ug MMC","WT E10 MMC","WT E10 MMC"),con_name="WT +MMC")

all_data[["dDBDdCTDnoMMC"]] <- merge_dataset(datasets,c("dDBDdCTD G9 5nM noMMC","dDBDcCTD G9 noMMC","dDBDdCTD G9 noMMC"),con_name="dDBDdCTD -MMC")

all_data[["dDBDdCTDMMC"]] <- merge_dataset(datasets,c("dDBDdCTD G9 5nM 1ug MMC","dDBDdCTD G9 MMC","dDBDdCTD G9 MMC"),con_name="dDBDdCTD +MMC")

all_data[["dDBDnoMMC"]] <- merge_dataset("D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD",c("dDBD A4 noMMC"),con_name="dDBD -MMC")

all_data[["dDBDMMC"]] <- merge_dataset("D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD",c("dDBD A4 MMC"),con_name="dDBD +MMC")

msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="HU",threshold=0.05,order=c(6,5,2,1,4,3))

###IR
all_data <- list()
datasets <- c("Y:/Maarten/Genetics/170517 BRCA2 tracking MMC","Y:/Maarten/Genetics/170623 BRCA2 MMC tracking","D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD")
all_data[["WTnoMMC"]] <- merge_dataset(datasets,c("WT E10 5nM noMMC","WT E10 noMMC2h","WT E10 noMMC"),con_name="WT -MMC")

all_data[["WTMMC"]] <- merge_dataset(datasets,c("WT E10 5nM 1ug MMC","WT E10 MMC","WT E10 MMC"),con_name="WT +MMC")

all_data[["dDBDdCTDnoMMC"]] <- merge_dataset(datasets,c("dDBDdCTD G9 5nM noMMC","dDBDcCTD G9 noMMC","dDBDdCTD G9 noMMC"),con_name="dDBDdCTD -MMC")

all_data[["dDBDdCTDMMC"]] <- merge_dataset(datasets,c("dDBDdCTD G9 5nM 1ug MMC","dDBDdCTD G9 MMC","dDBDdCTD G9 MMC"),con_name="dDBDdCTD +MMC")



##Old code
# directory <- "D:/Stack/Genetics/171208 N-BRCA2 Halo H8 MMC/"
# load(file.path(directory,"msd_fit_all.Rdata"))
# msd_fit_all$`-MMC`$cellID <- paste0('171208_',msd_fit_all$`-MMC`$cellID)
# msd_fit_all$`+MMC`$cellID <- paste0('171208_',msd_fit_all$`+MMC`$cellID)
#
# all$`-MMC` <- msd_fit_all$`-MMC`
# all$`+MMC` <- msd_fit_all$`+MMC`
#
# directory <- "D:/Stack/Genetics/171201 N-BRCA2 Halo H8 MMC/"
# load(file.path(directory,"msd_fit_all.Rdata"))
# msd_fit_all$`-MMC`$cellID <- paste0('171201_',msd_fit_all$`-MMC`$cellID)
# msd_fit_all$`+MMC`$cellID <- paste0('171201_',msd_fit_all$`+MMC`$cellID)
#
# all$`-MMC` <- rbind(all$`-MMC`,msd_fit_all$`-MMC`)
# all$`+MMC` <- rbind(all$`+MMC`,msd_fit_all$`+MMC`)
#
# directory <- "D:/Stack/Genetics/171117 N-BRCA2 MMC/"
# load(file.path(directory,"msd_fit_all.Rdata"))
# msd_fit_all$`H8 -MCC`$cellID <- paste0('171117_',msd_fit_all$`H8 -MCC`$cellID)
# msd_fit_all$`H8 +MCC`$cellID <- paste0('171117_',msd_fit_all$`H8 +MCC`$cellID)
#
# all$`-MMC` <- rbind(all$`-MMC`,msd_fit_all$`H8 -MCC`)
# all$`+MMC` <- rbind(all$`+MMC`,msd_fit_all$`H8 +MCC`)
#
# all$`-MMC`$condition <- "-MMC"
# all$`+MMC`$condition <- "+MMC"
# msd_fit_all <-all
