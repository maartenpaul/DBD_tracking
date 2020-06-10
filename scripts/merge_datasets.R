library(MSDtracking)

# BRCA2 HU ----------------------------------------------------------------
all_data <- list()
datasets <- c("D:/Stack/Genetics/180110 BRCA2 WTG10_dDdC_F4",
              "O:/Maarten/Genetics/170711 BRCA2 dDBD CTD Halo HU",
              "O:/Maarten/Genetics/170523 BRCA2 tracking HU",
              "D:/Stack/Genetics/180502 BRCA2 Halo dDBD+CTD HU/",
              "D:/Stack/Genetics/180606 BRCA2 WT dCTD Halo/")

all_data[["WT -HU"]] <- merge_dataset(datasets[c(2,3,4,5)],c("WT G10 noHU","WT E10 noHU","WT E10 noHU","WT G10 -HU","WT G10 -HU")[c(2,3,4,5)],con_name="WT -HU")

all_data[["WT +HU"]] <- merge_dataset(datasets[c(2,3,4,5)],c("WT G10 +HU 2h","WT E10 HU","WT E10 HU","WT G10 +HU","WT G10 +HU")[c(2,3,4,5)],con_name="WT +HU")

all_data[["dDBDdCTD -HU"]] <- merge_dataset(datasets[c(2,3,4)],c("dDBDdCTD F4 noHU","dDBDdCTD G9 noHU","dDBD G9 noHU","dDBDdCTD F4 -HU")[c(2,3,4)],con_name="dDBDdCTD -HU")

all_data[["dDBDdCTD +HU"]] <- merge_dataset(datasets[c(2,3,4)],c("dDBDdCTD F4 +HU 2h","dDBDdCTD G9 HU","dDBD G9 HU","dDBDdCTD F4 +HU")[c(2,3,4)],con_name="dDBDdCTD +HU")

all_data[["dDBD -HU"]] <- merge_dataset(datasets[2],c("dDBD+CTD A4 noHU"),con_name="dDBD -HU")

all_data[["dDBD +HU"]] <- merge_dataset(datasets[2],c("dDBD+CTD A4 HU"),con_name="dDBD +HU")

all_data[["dCTD -HU"]] <- merge_dataset(datasets[5],c("dCTD A2 -HU"),con_name="dCTD -HU")

all_data[["dCTD +HU"]] <- merge_dataset(datasets[5],c("dCTD A2 +HU"),con_name="dCTD +HU")

save(all_data,file = file.path("D:/Stack/Genetics/Tracking data merged/","BRCA2_HU_all_data.Rdata"))
msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_HU",threshold=0.05,merge_var = "experiment")
msd_histogram(all_data[c(1,2)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_HU_WT",threshold=0.05)
msd_histogram(all_data[c(3,4)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_HU_dDBDdCTD",threshold=0.05)
msd_histogram(all_data[c(5,6)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_HU_dDBD",threshold=0.05)
msd_histogram(all_data[c(1,2,3,4)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_HU+dDBDdCTD_WT",threshold=0.05)

all_data2 <- ldply(all_data)
cond <- data.frame(ldply(strsplit(all_data2$condition,split = " ")))
names(cond) <- c("type","treatment")
all_data2 <- cbind(all_data2,cond)

ggplot(all_data2,aes(x=D,group=treatment,col=treatment))+geom_histogram(binwidth=0.2,position = "identity",aes(y=..count../sum(..count..)),alpha=0.3)+scale_x_log10()+facet_grid(.~type)

ggplot(all_data2,aes(x=D,col=treatment))+scale_x_log10()+facet_grid(.~type)+geom_histogram(binwidth=0.2,aes(y=..density..), alpha=0.3,position="identity")

# BRCA2 MMC ---------------------------------------------------------------
all_data <- list()
datasets <- c("F://170517 BRCA2 tracking MMC",
              "F://170623 BRCA2 MMC tracking",
              "D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD")
all_data[["WT -MMC"]] <- merge_dataset(datasets[c(1,2,3)],c("WT E10 5nM noMMC","WT E10 noMMC2h","WT E10 noMMC")[c(1,2,3)],con_name="WT -MMC")

all_data[["WT +MMC"]] <- merge_dataset(datasets[c(1,2,3)],c("WT E10 5nM 1ug MMC","WT E10 MMC","WT E10 MMC")[c(1,2,3)],con_name="WT +MMC")

all_data[["dDBDdCTD -MMC"]] <- merge_dataset(datasets[c(1,2,3)],c("dDBDdCTD G9 5nM noMMC","dDBDcCTD G9 noMMC","dDBDdCTD G9 noMMC")[c(1,2,3)],con_name="dDBDdCTD -MMC")

all_data[["dDBDdCTD +MMC"]] <- merge_dataset(datasets[c(1,2,3)],c("dDBDdCTD G9 5nM 1ug MMC","dDBDdCTD G9 MMC","dDBDdCTD G9 MMC")[c(1,2,3)],con_name="dDBDdCTD +MMC")

all_data[["dDBD -MMC"]] <- merge_dataset("D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD",c("dDBD A4 noMMC"),con_name="dDBD -MMC")

all_data[["dDBD +MMC"]] <- merge_dataset("D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD",c("dDBD A4 MMC"),con_name="dDBD +MMC")

save(all_data,file = file.path("D:/Stack/Genetics/Tracking data merged/","BRCA2_MMC_all_data.Rdata"))
msd_histogram(all_data,directory = file.path("D:/Temp_data/"),name="MMC",threshold=0.05)
msd_histogram(all_data[c(1,2)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_MMC_WT",threshold=0.05)
msd_histogram(all_data[c(3,4)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_MMC_dDBDdCTD",threshold=0.05)
msd_histogram(all_data[c(1,2,3,4)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_MMC_WT+dDBDdCTD",threshold=0.05)
msd_histogram(all_data[c(5,6)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_MMC_dDBD",threshold=0.05)

# BRCA2 IR ----------------------------------------------------------------

all_data <- list()
datasets <- c("D:/Stack/Genetics/170510 particle tracking BRCA2/","F://180117 BRCA2 dDdC IR tracking")
all_data[["WT -IR"]] <- merge_dataset(datasets,c("BRCA2 WT Halo E10 noIR","WT G10 -IR"),con_name="WT -IR")

all_data[["WT +IR"]] <- merge_dataset(datasets,c("BRCA2 WT Halo E10 5gyIR 2h","WT G10 +IR"),con_name="WT +IR")

all_data[["dDBDdCTD -IR"]] <- merge_dataset(datasets,c("BRCA2 dDBDdCTD Halo G9 noIR","dDBDdCTD F4 -IR"),con_name="dDBDdCTD -IR")

all_data[["dDBDdCTD +IR"]] <- merge_dataset(datasets,c("BRCA2 dDBDdCTD Halo G9 5gyIR 2h","dDBDdCTD F4 +IR"),con_name="dDBDdCTD +IR")

save(all_data,file = file.path("D:/Stack/Genetics/Tracking data merged/","BRCA2_IR_all_data.Rdata"))
load(file.path("D:/Stack/Genetics/Tracking data merged/","BRCA2_IR_all_data.Rdata"))

msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_IR",threshold=0.05)
msd_histogram(all_data[c(1,2)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_WT_IR",threshold=0.05)
msd_histogram(all_data[c(3,4)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_dDBDdCTD_IR",threshold=0.05)

# BRCA2 IR iRFP720 ----------------------------------------------------------------

all_data <- list()
datasets <- c("F:/181005 BRCA2-Halo PCNA-iRFP IR tracking/","F:/181010 BRCA2-Halo PCNA-iRFP IR tracking/","F:/181012 BRCA2-Halo PCNA-iRFP IR tracking/")
all_data[["WT -IR"]] <- merge_dataset(datasets[c(1,3)],c("WT G10 -IR","WT G10 -IR"),con_name="WT -IR")

all_data[["WT +IR"]] <- merge_dataset(datasets[c(1,3)],c("WT G10 +IR","WT G10 +IR"),con_name="WT +IR")

all_data[["dDBD -IR"]] <- merge_dataset(datasets[c(1,2)],c("dDBD E4 -IR","dDBD E4 -IR"),con_name="dDBD -IR")

all_data[["dDBD +IR"]] <- merge_dataset(datasets[c(1,2)],c("dDBD E4 +IR","dDBD E4 +IR"),con_name="dDBD +IR")

all_data[["dCTD -IR"]] <- merge_dataset(datasets[c(1,2)],c("dCTD A2 -IR","dCTD A2 -IR"),con_name="dCTD -IR")

all_data[["dCTD +IR"]] <- merge_dataset(datasets[c(1,2)],c("dCTD A2 +IR","dCTD A2 +IR"),con_name="dCTD +IR")

all_data[["dDBDdCTD -IR"]] <- merge_dataset(datasets[c(2,3)],c("dDBDdCTD F4 -IR","dDBDdCTD F4 -IR"),con_name="dDBDdCTD -IR")

all_data[["dDBDdCTD +IR"]] <- merge_dataset(datasets[c(2,3)],c("dDBDdCTD F4 +IR","dDBDdCTD F4 +IR"),con_name="dDBDdCTD +IR")


save(all_data,file = file.path("D:/Stack/Genetics/Tracking data merged/","BRCA2_IR_all_data.Rdata"))
load(file.path("D:/Stack/Genetics/Tracking data merged/","BRCA2_IR_all_data.Rdata"))

msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_IR",threshold=0.05)#,merge_var = "cellID")
msd_histogram(all_data[c(7,8)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_WT_IR",threshold=0.1)#,merge_var = "cellID")
msd_histogram(all_data[c(3,4)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_dDBDdCTD_IR",threshold=0.05,merge_var = "cellID")

# BRCA2 dDdC#2 IR iRFP720 ----------------------------------------------------------------

all_data <- list()
datasets <- c("D:/Imaging data/190206 MPexp1902_01 dDdC IR tracking/","D:/Imaging data/190206_2 MPexp1902_01 dDdC IR tracking_2/")
all_data[["dDBDdCTD -IR"]] <- merge_dataset(datasets,c("dDdC F4 -IR","dDdC F4 -IR"),con_name="dDBDdCTD -IR")

all_data[["dDBDdCTD +IR"]] <- merge_dataset(datasets,c("dDdC F4 +IR","dDdC F4 +IR"),con_name="dDBDdCTD +IR")


save(all_data,file = file.path("D:/Stack/Genetics/Tracking data merged/","BRCA2_dDdCv2_IR_all_data.Rdata"))
load(file.path("D:/Stack/Genetics/Tracking data merged/","BRCA2_dDdCv2_IR_all_data.Rdata"))

msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_dDdCv2_IR",threshold=0.05)#,merge_var = "cellID")



# BRCA2 N-terminal iRFP720  --------------------------------------------------------------

all_data <- list()
datasets <- c("F:/181210 Halo-BRCA2 iRFP270 tracking/","F:/181212 Halo-BRCA2 iRFP270 tracking/","F:/181214 PALB2-Halo iRFP270-PCNA tracking")
all_data[["WT -IR"]] <- merge_dataset(datasets[2],c("Halo-BRCA2 F4 -IR"),con_name="WT -IR")

all_data[["WT +IR"]] <- merge_dataset(datasets[2],c("Halo-BRCA2 F4 +IR"),con_name="WT +IR")

all_data[["d1-40 -IR"]] <- merge_dataset(datasets[2],c("d1-40 E3 -IR"),con_name="d1-40 -IR")

all_data[["d1-40 +IR"]] <- merge_dataset(datasets[2],c("d1-40 E3 +IR"),con_name="d1-40 +IR")

all_data[["PALB2 -IR"]] <- merge_dataset(datasets[3],c("PALB2-Halo A6 -IR"),con_name="PALB2 -IR")

all_data[["PALB2 +IR"]] <- merge_dataset(datasets[3],c("PALB2-Halo A6 +IR"),con_name="PALB2 +IR")

save(all_data,file = file.path("D:/Stack/Genetics/Tracking data merged/","N_BRCA2_IR_iRFP_all_data.Rdata"))
load(file.path("D:/Stack/Genetics/Tracking data merged/","N_BRCA2_IR_iRFP_all_data.Rdata"))

msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="N_BRCA2_iRFP",threshold=0.05)#,merge_var = "experiment")
#msd_histogram(all_data[c(1,2)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_WT_IR",threshold=0.1)#,merge_var = "cellID")
#msd_histogram(all_data[c(3,4)],directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRCA2_dDBDdCTD_IR",threshold=0.05,merge_var = "cellID")




# BRCA2 MMC iRFP720-PCNA --------------------------------------------------


all_data <- list()
# 180702 WT (only 5 +MMC) dCTD D:/OneDrive/Data2/180702 BRCA2_Halo_WT_dCTD PCNA MMC
#
# Exp1905_01 WT and dDBDdCTD (~11 cells) D:/OneDrive/Data2/190504 MPexp1905_01 BRCA2 tracking MMC
#
# Exp1912_003 Run 1 dDBD and dDBDdCTD (~9 cells) D:/OneDrive/Data2/191211 ExpMP1912_004 MMC tracking/Run 1
#
# Exp1912_003 Run 2 dDBD and dCTD (~9 cells) D:/OneDrive/Data2/191211 ExpMP1912_004 MMC tracking/Run 2


datasets <- c("D:/OneDrive/Data2/180702 BRCA2_Halo_WT_dCTD PCNA MMC/","D:/OneDrive/Data2/190504 MPexp1905_01 BRCA2 tracking MMC/",
              "D:/OneDrive/Data2/191211 ExpMP1912_004 MMC tracking/Run 1/","D:/OneDrive/Data2/191211 ExpMP1912_004 MMC tracking/Run 2/","D:/OneDrive/Data2/200207 ExpMP2002_002 MMC tracking/")
all_data[["WT -MMC"]] <- merge_dataset(datasets[c(1,2,5)],c("WT G10 mock","WTG10 iRFP720PCNA mock","WT G10 -MMC"),con_name="WT -MMC")

all_data[["WT +MMC"]] <- merge_dataset(datasets[c(1,2,5)],c("WT G10 MMC","WTG10 iRFP720PCNA MMC","WT G10 -MMC"),con_name="WT +MMC")

all_data[["dDBD -MMC"]] <- merge_dataset(datasets[c(3,4)],c("dDBD E4 -MMC","dDBD E4 -MMC"),con_name="dDBD -MMC")

all_data[["dDBD +MMC"]] <- merge_dataset(datasets[c(3,4)],c("dDBD E4 +MMC","dDBD E4 +MMC"),con_name="dDBD +MMC")

all_data[["dCTD -MMC"]] <- merge_dataset(datasets[c(1,4)],c("dCTD A2 mock","dCTD A2 -MMC"),con_name="dCTD -MMC")

all_data[["dCTD +MMC"]] <- merge_dataset(datasets[c(1,4)],c("dCTD A2 MMC","dCTD A2 +MMC"),con_name="dCTD +MMC")

all_data[["dDBDdCTD -MMC"]] <- merge_dataset(datasets[c(2,3,5)],c("dDBDdCTD F4 iRFP720PCNA mock","dDBDdCTD F4 -MMC","dDdC F4 -MMC"),con_name="dDBDdCTD -MMC")

all_data[["dDBDdCTD +MMC"]] <- merge_dataset(datasets[c(2,3,5)],c("dDBDdCTD F4 iRFP720PCNA MMC","dDBDdCTD F4 +MMC","dDdC F4 +MMC"),con_name="dDBDdCTD +MMC")


save(all_data,file = file.path("D:/OneDrive/Data2/Tracking data merged/","BRCA2_MMC_PCNA_all_data.Rdata"))
load(file.path("D:/OneDrive/Data2/Tracking data merged/","BRCA2_MMC_PCNA_all_data.Rdata"))

msd_histogram(all_data,directory = file.path("D:/OneDrive/Data2/Tracking data merged/"),name="BRCA2_MMC_PCNA",threshold=0.05)#,merge_var = "cellID")
msd_histogram(all_data[c(7,8)],directory = file.path("D:/OneDrive/Data2/Tracking data merged/"),name="BRCA2_WT_MMC_PCNA",threshold=0.1)#,merge_var = "cellID")
msd_histogram(all_data[c(3,4)],directory = file.path("D:/OneDrive/Data2/Tracking data merged/"),name="BRCA2_dDBDdCTD_IR",threshold=0.05,merge_var = "cellID")

# PALB2 iRFP720  --------------------------------------------------------------
all_data <- list()
datasets <- c("F:/181214 PALB2-Halo iRFP270-PCNA tracking")
all_data[["PALB2 -IR"]] <- merge_dataset(datasets[1],c("PALB2-Halo A6 -IR"),con_name="PALB2 -IR")

all_data[["PALB2 +IR"]] <- merge_dataset(datasets[1],c("PALB2-Halo A6 +IR"),con_name="PALB2 +IR")


save(all_data,file = file.path("D:/Stack/Genetics/Tracking data merged/","PALB2_IR_iRFP_all_data.Rdata"))
load(file.path("D:/Stack/Genetics/Tracking data merged/","PALB2_IR_iRFP_all_data.Rdata"))

msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="PALB2_IR_iRFP",threshold=0.05)


# ##PALB2 HU --------------------------------------------------------------
all_data <- list()
stop("TO DO")
datasets <- c("F://170517 BRCA2 tracking MMC","F://170623 BRCA2 MMC tracking","D:/Stack/Genetics/170809 tracking MMC WT_dDBD_+CTD")
all_data[["WTnoMMC"]] <- merge_dataset(datasets,c("WT E10 5nM noMMC","WT E10 noMMC2h","WT E10 noMMC"),con_name="WT -MMC")

all_data[["WTMMC"]] <- merge_dataset(datasets,c("WT E10 5nM 1ug MMC","WT E10 MMC","WT E10 MMC"),con_name="WT +MMC")


# PALB2 MMC ---------------------------------------------------------------
all_data <- list()
datasets <- c("D:/Stack/Genetics/171024 PALB2-Halo tracking")
all_data[["PALB2 -MMC"]] <- merge_dataset(datasets,c("PALB2 F6 mock"),con_name="PALB2 -MMC")

all_data[["PALB2 +MMC"]] <- merge_dataset(datasets,c("PALB2 F6 MMC 2h"),con_name="PALB2 +MMC")

msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="PALB2 MMC",threshold=0.05)


# Halo-BRCA2 WT d1-40 -----------------------------------------------------
all_data <- list()
datasets <- c("D:/Stack/Genetics/171130 d1-40 Halo-BRCA2","D:/Stack/Genetics/171130 N-BRCA2 Halo H8 HU")
all_data[["Halo-BRCA2"]] <- merge_dataset(datasets[2],c("-HU"),con_name="Halo-BRCA2")

all_data[["Halo-BRCA2 d1-40"]] <- merge_dataset(datasets[1],c("F9"),con_name="Halo-BRCA2 d1-40")

msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="HaloTag BRCA2 d1-40",threshold=0.05)




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



# H2B/NLS ---------------------------------------------------------------------
all_data <- list()
#all_data[["H2B"]] <- merge_dataset(datasets[1],c("H2B"),con_name="H2B")
load("Y:/Maarten/Genetics/170809 tracking MMC WT_dDBD_+CTD/H2B/msd_fit.Rdata")
all_data[["H2B"]] <- ldply(tracks)
all_data[["H2B"]]$condition <- "H2B"

load("Y:/Maarten/Genetics/170512 BRCA2 slow tracking/Halo-NLS/30 ms/msd_fit.Rdata")
all_data[["NLS"]] <- ldply(tracks)
all_data[["NLS"]]$condition <- "NLS"
msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="H2B-NLS",threshold=0.05)


#BRC3 preliminary
all_data <- list()
datasets <- c("D:/Stack/Genetics/180502 BRCA2 Halo dDBD+CTD HU","D:/Stack/Genetics/180502 BRCA2 Halo WT BRC3-F2A-GFP",
              "D:/Stack/Genetics/180426 BRCA2 Halo dDBD+CTD HU","D:/Stack/Genetics/180426 BRC3-F2A-GFP/WT G10 stable #2")
all_data[["WT -BRC3"]] <- merge_dataset(datasets[c(1,3)],c("WT G10 -HU","BRCA2 WT G10 -HU"),con_name="WT -BRC3")

all_data[["WT +BRC3"]] <- merge_dataset(datasets[c(2,4)],c("+BRC3","+BRC3"),con_name="WT +BRC3")
msd_histogram(all_data,directory = file.path("D:/Stack/Genetics/Tracking data merged/"),name="BRC3",threshold=0.05)

