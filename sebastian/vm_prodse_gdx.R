#script reads a bunch of fulldata.gdxs and writes the vm_prodse 

library(gdx)

wd_dir <- "C:/Users/rauner/Documents/PIK/fasstr_sebastian"
setwd(wd_dir)
dir_root           <- paste0(getwd()    ,'/mnt/FASST_IDL_PACK/')              # main folder
remind_output      <- paste0(getwd()    ,'/mnt/FASST_IDL_PACK/REMIND_OUTPUT/internalization')
dir_outtab         <- paste0(dir_root,  'OUTPUT/TABLES')                      # output (tables)


scenario_folders <- paste0(list.dirs(paste0(remind_output), recursive = FALSE),'/')

#initialize
out <- read.gdx(paste0(scenario_folders[1], 'fulldata.gdx'),requestList.name='vm_prodSe', fields= 'l', squeeze = F)

for( scenario_folder in scenario_folders) {
  
  #for the scenario name
  remind_config      <- paste0(scenario_folder,     'config.Rdata')
  load(remind_config)
  scenario <- sapply(strsplit(cfg$gms$c_expname,'-'), "[[", 2)
  
  AP_scenario <- cfg$gms$cm_APscen
  
  if (scenario_folder == scenario_folders[1]) {
    out <- read.gdx(paste0(scenario_folders[1], 'fulldata.gdx'),requestList.name='vm_prodSe', fields= 'l', squeeze = F)
    out$'scenario' <- scenario
  } else {
    se <- read.gdx(paste0(scenario_folder, 'fulldata.gdx'),requestList.name='vm_prodSe', fields= 'l', squeeze = F)
    se$'scenario' <- scenario
    out <- rbind(out,se)
  }
  
}

#rearange columns

out <- out %>%
  select(scenario, everything()) %>%
  select(-"all_enty",-"all_enty.1")


# if zero put small number
out$value[out$value==0] <- 10^-9

write.csv(out,file=paste(dir_outtab,'/vm_prodse_internalization.csv',sep=''))
