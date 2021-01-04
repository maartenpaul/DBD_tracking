library(MSDtracking)

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
