#calculate the direct economic cost

library(dplyr)
library(tidyr)

#setwd("C:/Users/rauner/Documents/PIK/fasstr_sebastian")

wd_dir <- "C:/Users/rauner/Documents/PIK/fasstr_sebastian"
setwd(wd_dir)

dir_root   <- paste0(getwd(),   '/mnt/FASST_IDL_PACK/')            # main folder
dir_ancil  <- paste0(dir_root,  'ANCILLARY')                       # ancillary data
dir_DALYS  <- paste0(dir_ancil, '/DALYS')                          # DALYS
dir_outtab <- paste0(dir_root,  'OUTPUT/TABLES')           # output (tables)

# Direct economic cost of air pollution health effects
# As a cost per Daly
# Productivity loss due to loss of working days
# from “Holland, M. (2014), Cost-benefit Analysis of Final Policy Scenarios for the EU Clean Air Package,
# http://ec.europa.eu/environment/air/pdf/TSAP%20CBA.pdf “
u_direct        <- 256127 / 25



# 
# #read the premature deaths to DALY ratio
# daly_ratio           <- read.csv(file=file.path(dir_DALYS, 'daly_ratio.csv'),sep=';')
# 
# daly_ratio[is.na(daly_ratio)] <- 0
# 
# #read the results table
# write_out_data_frame <- read.csv(file.path(dir_outtab,paste0('results_',RUNSUFFIX, '.csv')),sep=',')
# 
# 
# #calcualte the DALYS and write them in the results table
# #alri are <5, the rest >30
# 
# #join results and DALYS
# results <- inner_join(write_out_data_frame,daly_ratio, by=c("iso.iso.y"="ISO")) %>% mutate(daly_copd = data_df_mort_pm.copd1* copd) %>%
#                                                                          mutate(daly_alri = data_df_mort_pm.alri1* alri) %>%
#                                                                          mutate(daly_ihd  = data_df_mort_pm.ihd1 * ihd)  %>%
#                                                                          mutate(daly_stroke = data_df_mort_pm.stroke1* stroke)%>%
#                                                                          mutate(daly_lc = data_df_mort_pm.lc1* lc)%>%
#                                                                          mutate(daly_total = daly_copd+daly_alri+
#                                                                                              daly_ihd +daly_lc+
#                                                                                              daly_stroke)

#calculate the new social cost of DALYS
source(paste0(getwd(),   '/functions/calcMonetization.R')   )
# calculate the DALY coefficients
DALY = u_direct * calcMonetization()
DALY[is.na(DALY)] <- 0
DALY <- as.data.frame(DALY) %>% dplyr::select('Region', 'Year', 'Value')
names(DALY)[3] <- 'direct_cost_per_daly'

DALY$Year <- as.character(DALY$Year)

DALY$Year <- as.integer(DALY$Year)


#multiply dalys with results dalys
results <- inner_join(results, DALY, by=c("iso.iso.y"="Region", 'year'='Year')) %>% mutate(direct_daly_cost = daly_total * direct_cost_per_daly)
results$'scenario' <- scenario
results[,1:2]<- list(NULL)
results$daly_cost <- results$indirect_daly_cost + results$direct_daly_cost

write.table(results, file.path(dir_outtab,paste0('results_dalys_',RUNSUFFIX, '.csv')),sep=';', row.names = F)


