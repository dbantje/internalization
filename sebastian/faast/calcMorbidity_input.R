#calculate the YLL and DALY per country and writes them in a csv which is used in the DALY part
# read from GBD source
# GBD Results tool:
#   Use the following to cite data included in this download:
#   Global Burden of Disease Collaborative Network.
# Global Burden of Disease Study 2016 (GBD 2016) Results.
# Seattle, United States: Institute for Health Metrics and Evaluation (IHME), 2017.
# Available from http://ghdx.healthdata.org/gbd-results-tool.

# extract DALY/deaths for the 5 endpoints for <5 >30 age groupe and iso country

library(dplyr)
library(tidyr)

#setwd("C:/Users/rauner/Documents/PIK/fasstr_sebastian")

wd_dir <- "C:/Users/rauner/Documents/PIK/fasstr_sebastian"
setwd(wd_dir)

dir_root   <- paste0(getwd(),   '/mnt/FASST_IDL_PACK/')            # main folder
dir_ancil  <- paste0(dir_root,  'ANCILLARY')                       # ancillary data
dir_DALYS  <- paste0(dir_ancil, '/DALYS')                          # DALYS

# 
# #read location id
# location_GBD <- read.table(file=file.path(dir_DALYS, "IHME_GBD_2016_ALL_LOCATION_HIERARCHIES_Y2018M07D20.csv"), sep = ';', header = TRUE)
# location_GBD <- location_GBD %>% filter(level == 2)
# 
# # delte duplication
# 
# 
# 
# #read REMIND mapping
# regionmapping_remind <- read.table(file=file.path(dir_DALYS, "regionmappingGBD.csv"), sep = ';', header = TRUE)
# 
# 
# spatial_mapping <- stringdist_left_join(regionmapping_remind, location_GBD, by=c('X'='location_name'), method = 'osa')
# 

#read regionmappingGBD
regionmappingGBD <-  read.table(file=file.path(dir_DALYS, "regionmappingGBD.csv"), sep = ';', header = TRUE)


#read source
input  <- rbind(read.table(file=file.path(dir_DALYS, "IHME-GBD_2017_DATA-14a3b68e-1.csv"), sep = ',', header = T),
          read.table(file=file.path(dir_DALYS, "IHME-GBD_2017_DATA-14a3b68e-2.csv"), sep = ',', header = T))




# filter for DALYS and deaths, number and both


deaths <- input %>% filter(measure_id == 1) %>% filter(metric_id == 1) %>% filter(sex_id == 3)
dalys <-  input %>% filter(measure_id == 2) %>% filter(metric_id == 1) %>% filter(sex_id == 3)

daly_ratio <- inner_join(dalys,deaths, by=c("location_id", "age_id", "cause_id"))

#calculate the ratio between both
daly_ratio$daly_ratio <- daly_ratio$val.x/daly_ratio$val.y

#select relevant data
daly_ratio <- daly_ratio %>% select(measure_id.x, location_id,sex_id.x, age_name.x , cause_name.x, metric_id.x, daly_ratio)




#from long to wide age_id

daly_ratio <- spread(daly_ratio, age_name.x, daly_ratio)

#aggregate to >30
#add new lines with >30
#list of >30 ages
#22 minus 1, 6, 7, 8, 9, 10 
daly_ratio %>%
dplyr::mutate('>30' = 1+22)

#that is not really good here

daly_ratio$'>30' <- (daly_ratio$'30 to 34'+
                       daly_ratio$'35 to 39'+
                       daly_ratio$'40 to 44'+
                       daly_ratio$'45 to 49'+
                       daly_ratio$'50 to 54'+
                       daly_ratio$'55 to 59'+
                       daly_ratio$'60 to 64'+
                       daly_ratio$'65 to 69'+
                       daly_ratio$'70 to 74'+
                       daly_ratio$'75 to 79'+
                       daly_ratio$'80 to 84'+
                       daly_ratio$'85 to 89'+
                       daly_ratio$'95 plus')/13

#aggregate to ISO countries

daly_ratio$ISO <- regionmappingGBD$CountryCode[plyr::mapvalues(daly_ratio$location_id, from = regionmappingGBD$location_id, to = regionmappingGBD$ISO)]

#select relevant cols
daly_ratio <- daly_ratio %>% select("ISO","cause_name.x", "Under 5", ">30")

#rename relevant cols
names(daly_ratio)[3]<-"frac_05"
names(daly_ratio)[4]<-"frac_30"

#rename cause_name to the fasst output names
daly_ratio$cause_name.x <- as.character(daly_ratio$cause_name.x)
daly_ratio$cause_name.x[daly_ratio$cause_name.x  == "Ischemic heart disease"] <- "ihd"
daly_ratio$cause_name.x[daly_ratio$cause_name.x  == "Chronic obstructive pulmonary disease"] <- "copd"
daly_ratio$cause_name.x[daly_ratio$cause_name.x  == "Stroke"] <- "stroke"
daly_ratio$cause_name.x[daly_ratio$cause_name.x  == "Tracheal, bronchus, and lung cancer"] <- "lc"
daly_ratio$cause_name.x[daly_ratio$cause_name.x  == "Lower respiratory infections"] <- "alri"
daly_ratio$cause_name.x[daly_ratio$cause_name.x  == "All causes"] <- "all_mort"


#long to wide
daly_ratio <- reshape(daly_ratio, idvar = "ISO", timevar = "cause_name.x", direction = "wide")

#only select relevant columns
daly_ratio <- daly_ratio %>% select('ISO','frac_05.alri','frac_30.ihd','frac_30.copd' , 'frac_30.stroke','frac_30.lc' )

#rename columns
names(daly_ratio)[2]<-"alri"
names(daly_ratio)[3]<-"ihd"
names(daly_ratio)[4]<-"copd"
names(daly_ratio)[5]<-"stroke"
names(daly_ratio)[6]<-"lc"


write.csv(daly_ratio,file=file.path(dir_DALYS, 'daly_ratio.csv'))







