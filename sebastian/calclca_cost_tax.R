#disaggregate the LCA to iso country and get dalys and cost from AP table results

# disaggregate the ap dalys to sectors with LCA ap ratio

#read lca results, read talbe ap results

# calculates the cost per technology and region and tax

library(dplyr)
library(tidyr)
library(madrat)
library(mrcommons)
library(ggplot2)
library(quitte)

wd_dir <- "C:/Users/rauner/Documents/PIK/fasstr_sebastian"
setwd(wd_dir)

do.plots <- F

dir_root           <- paste0(getwd()    ,'/mnt/FASST_IDL_PACK/')              # main folder
dir_ancil  <- paste0(dir_root,  'ANCILLARY')                       # ancillary data
dir_DALYS  <- paste0(dir_ancil, '/DALYS')                          # DALYS
dir_outtab <- paste0(dir_root,  'OUTPUT/TABLES/')           # output (tables)
#dir_outtab <- paste0(dir_root,  'OUTPUT/TABLES/enavi/')           # output (tables)
dir_cost   <- paste0(dir_root,  'OUTPUT/cost/')           # output (tables)
dir_lca    <- paste0('C:/Users/rauner/Documents/PIK/holi_coal_exit/rmnd-lca/script/results/')            # main folder
dir_plot_out       <- paste0(getwd(),  '/plots/')                              # output
dir_SE    <- paste0('C:/Users/rauner/Documents/PIK/holi_coal_exit/rmnd-lca/script/data/Remind output files/')   


#scenarios
scenario_files <- list.files(path = dir_outtab, pattern = 'results_dalys_SSP2')
lca_files      <- list.files(path = dir_lca, pattern = 'lca_REMIND_monetization_project.csv')

end_category <- c('ED_H', 'HH_H', 'RA_H',
                  'Policy Cost|Consumption Loss (billion US$2005/yr)',
                  'Policy Cost|GDP Loss (billion US$2005/yr)',
                  'Policy Cost|Additional Total Energy System Cost (billion US$2005/yr)',
                  'Policy Cost|Consumption + Current Account Loss (billion US$2005/yr)')
end_category_clear <- c('Ecosystem Damages [US$2005 from species.yr]',
                        'Human Health [US$2005 from DALY]',
                        'Ressource Depletion [US$2005]',
                        'Consumption Loss [US$2005/yr]',
                        'GDP Loss [US$2005/yr]',
                        'Additional Total Energy System Cost [US$2005/yr]',
                        'Consumption + Current Account Loss [US$2005/yr]')


category <- c(
              'PMF',
              'CC',
              'OD',
              'TA',
              'FE',
              'HT',
              'POF',
              
              'TET',
              'FET',
              'MET',
              'IR',
              'ALO',
              'ULO',
              'NLT',
              'FD',
              'MD',
              'ME',
              'WD',
              'Policy Cost|Consumption Loss (billion US$2005/yr)',
              'Policy Cost|GDP Loss (billion US$2005/yr)',
              'Policy Cost|Additional Total Energy System Cost (billion US$2005/yr)',
              'Policy Cost|Consumption + Current Account Loss (billion US$2005/yr)')

category_clear <- c('Air Pollution',
                    'Climate change',
                    'Ozone depletion',
                    'Terrestrial acidification',
                    'Freshwater eutrophication',
                    'Human toxicity',
                    'Air Pollution',
                    
                    'Terrestrial ecotoxicity',
                    'Freshwater ecotoxicity',
                    'Marine ecotoxicity',
                    'Ionising radiation',
                    'Agricultural land occupation',
                    'Urban land occupation',
                    'Natural land transformation',
                    'Fossil depletion',
                    'Metal depletion',
                    'Marine eutrophication',
                    'Water depletion',
                    'Consumption Loss',
                    'GDP Loss',
                    'Additional Total Energy System Cost',
                    'Consumption + Current Account Loss')

category_clear_impact <- c('Air Pollution',
                           'Climate change',
                           'Ozone depletion',
                           'Terrestrial acidification',
                           'Freshwater eutrophication',
                           'Human toxicity',
                           'Photochemical oxidant formation',
                           
                           'Terrestrial ecotoxicity',
                           'Freshwater ecotoxicity',
                           'Marine ecotoxicity',
                           'Ionising radiation',
                           'Agricultural land occupation',
                           'Urban land occupation',
                           'Natural land transformation',
                           'Fossil depletion',
                           'Metal depletion',
                           'Marine eutrophication',
                           'Water depletion')
# read and rbind tab results
ap <- do.call(rbind,
        lapply(paste0(dir_outtab,scenario_files), read.csv, sep=';'))

#read the lca results
lca <- read.csv(paste0(dir_lca,lca_files[1]), sep=';')
#lca <- read.csv(paste0(dir_lca,lca_files[2]), sep=';')


# filter out the totals and condense category and end category
lca <- lca %>% filter(category != 'total')
lca <- lca %>% filter(category != 'Unnamed: 0')


#take the 2015 values from Base for all other scenarios
lca_2015 <- lca %>% filter(scenario == 'Base') %>% filter(period == 2015)

lca  <- lca %>% filter(period != 2015)

for(scen in unique(lca$scenario)){
  
  lca_2015$scenario <- as.character(scen)
  lca <- rbind(lca,lca_2015)
}


# read the gdx SE

SE_gdx <- read.csv(paste0(dir_SE,'vm_prodse_internalization.csv'))
FE_gdx <- read.csv(paste0(dir_SE,'vm_prodfe_internalization.csv'))

SE_tech <- unique(SE_gdx$all_te)
FE_tech <- unique(FE_gdx$all_te)

names(FE_gdx)[names(FE_gdx) == 'ttot'] <- 'tall'

SE_gdx <- rbind(SE_gdx, FE_gdx)


#take the 2015 values from Base for all other scenarios
SE_gdx_2015 <- SE_gdx %>% filter(scenario == 'Base') %>% filter(tall == 2015)

SE_gdx  <- SE_gdx %>% filter(tall != 2015)

for(scen in unique(lca$scenario)){
  
  SE_gdx_2015$scenario <- as.character(scen)
  SE_gdx <- rbind(SE_gdx,SE_gdx_2015)
}



#add a total cost column
lca$cost_total_vsl <- lca$ED_H_cost_total_vsl + lca$RA_H_cost_total_vsl  + lca$HH_H_cost_total_vsl

# filter out world
lca <- lca %>% filter(region != 'World')

lca$region <- as.character(lca$region)

lca <-lca %>% filter(region != '')

#add here the policy cost and include them as an end_category in the lca

#filter the LCA impact

#SE <- lca  %>% select(c('scenario','period','region','category', 'end_category','tech','SE')) 
lca <- lca %>% select(c('scenario','period','region','category', 'end_category','tech','cost_total_vsl')) 
#year
year <- unique(lca$period)
#year <- 2015
lca <- lca %>% filter(period %in% year)
#SE  <- SE %>% filter(period %in% year)

names(lca)[names(lca) == "cost_total_vsl"] <- "value"
#names(SE)[names(SE) == "SE"] <- "value"

#xxx do we loose the ap?
# yes we do, fix it
# are we loosing duplicates here?
# urban land zero
# check fixing mid to end

# why is it so low, are we loosing duplicates here????


lca <- as.magpie(as.quitte(lca))
#SE <- as.magpie(as.quitte(SE))



#region mapping
region_mapping <- read.csv('regionmappingH12.csv', sep =';')
#region_mapping <- read.csv('regionmapping_22_EU11.csv', sep =';')

#pop as weights
pop                           <- calcOutput("Population",     years=year, aggregate = F ) 
pop=pop[,,"SSP2",pmatch=TRUE]
pop=pop[,2100,,invert=TRUE]
lca = madrat:::toolAggregate(lca, rel =  region_mapping, weight = pop)
#SE = speed_aggregate(SE, rel =  region_mapping, weight = pop)
lca = quitte::as.quitte(lca)
#SE = quitte::as.quitte(SE)
mitigation = read.csv(paste0(dir_cost,'/policy_cost_iso_trend.csv'), sep=',')
#mitigation = read.csv(paste0(dir_cost,'/policy_cost_iso_enavi.csv'), sep=',')

#select here only one of the mitigation measurement
#mitigation <- mitigation %>% filter(mitigation_variable == 'Policy Cost|Additional Total Energy System Cost (billion US$2005/yr)')
mitigation$mitigation_value <-mitigation$mitigation_value*10^9
#add columns to mitigation to be able to rbind it to LCA as a new end category
mitigation$model <- NA
mitigation$variable <- NA
mitigation$unit <- NA
mitigation$category <- mitigation$mitigation_variable
mitigation$end_category <- mitigation$mitigation_variable
mitigation$tech <- 'mitigation'
mitigation$value <- mitigation$mitigation_value

mitigation <- mitigation %>% select(paste(colnames(lca)))

lca <- rbind(lca,mitigation)

# does that really make sense here? its more a loss of ssectoral information.

#filter the AP endpoint of the LCA and make a relation of the techs out of it
# distribute the AP total dalys to the techs with this relation
lca_pmf <- lca %>% filter(category == 'PMF')
lca_o3 <- lca %>% filter(category == 'POF')
lca_pmf$value[is.na(lca_pmf$value)] <- 0
lca_o3$value[is.na(lca_o3$value)] <- 0

#add a total row
lca_pmf_total= lca_pmf %>% 
               dplyr::group_by(scenario, region, variable, unit,  period) %>% 
               dplyr::summarise(total = sum(value))


lca_o3_total= lca_o3 %>% 
  group_by(scenario, region, variable, unit,  period) %>% 
  dplyr::summarise(total = sum(value))

# join lca_pmf and total
lca_pmf <- dplyr::left_join(lca_pmf, lca_pmf_total)
lca_pmf$value <- lca_pmf$value / lca_pmf$total

lca_o3 <- dplyr::left_join(lca_o3, lca_o3_total)
lca_o3$value <- lca_o3$value / lca_o3$total



#filter only for dalys cost column of ap
ap_dalys <- ap %>% select('scenario','year','iso.iso.y', 'daly_cost')

# remove -rem-5 if its there

ap_dalys$scenario <- str_replace(ap_dalys$scenario, '-rem-5','')

#map ap to the lca region mapping



#rather distribute the ap with emissions as proxies?

#match lca_pmf and ap scenario names
#multiply the ap value with the lca_pmf value(ratio sector to total)
# there should be a bug

#remove the SSP2 from scenario
ap_dalys$scenario <- gsub("SSP2-",'',ap_dalys$scenario)


lca_pmf <- dplyr::left_join(lca_pmf, ap_dalys, by =c('scenario'='scenario',
                                               'period'  ='year',
                                               'region'  ='iso.iso.y'))
lca_pmf$value = lca_pmf$value * lca_pmf$daly_cost
lca_pmf <- lca_pmf %>% select(-daly_cost, -total)

lca_o3 <- dplyr::left_join(lca_o3, ap_dalys, by =c('scenario'='scenario',
                                                     'period'  ='year',
                                                     'region'  ='iso.iso.y'))
lca_o3$value = lca_o3$value * lca_o3$daly_cost
lca_o3 <- lca_o3 %>% select(-daly_cost, -total)

# we distribute the air pollution dalys with the PMF and POF ratios
# do a better job to distribute the AP stuff

# probably the problem is that Alois adjustes the fossils to lower their emission factors considerably
#but not the renewables

#1. ask Alois how he adjustes the air pollution emissions

#1. check with the activity manager where
#1. cut out the adjustement and rerun?
#1. adjust the photochemical ozone back to the old one



#2. do a better allocation, maybe not consider embodied emissions?
# distribute with emission

#0.5 for inequality

#lca_o3$value <- lca_o3$value *0.2
lca_o3$value <- lca_o3$value *0.0
lca_o3$category <- 'POF'

#
#lca_pmf$value <- lca_pmf$value *0.8
lca_pmf$value <- lca_pmf$value *0.0


# add this to the lca insted of particulate matter formation
# filter lca for non pmf and rbind lca and lca_pmf (adjusted with our ap results)

lca <- lca %>% filter(category != 'PMF' )
lca <- lca %>% filter(category != 'POF' )

lca <- rbind(lca, lca_pmf,lca_o3)

#make some REMIND region plots out of that
lca_remind <- dplyr::left_join(lca, region_mapping, by =c('region'='CountryCode'))
#SE <- dplyr::left_join(SE, region_mapping, by =c('region'='CountryCode'))

# add end_category and end_category_tech totals
# add category totals

#this total does not include the AP
lca_remind <- lca_remind %>% filter(category != 'total')

lca_remind <- lca_remind %>% select(-variable, -unit, -model)
lca_remind[lca_remind == NA] <- 0

#we dont include CC
lca_remind <- lca_remind %>% filter(category != 'CC')

lca_remind <- lca_remind %>% filter(category != 'FD')
lca_remind <- lca_remind %>% filter(category != 'MD')
# 
# #here we make the cost differential to the Base, because the policy cost are also relative to the Reference
# 
# # filter for Reference and the rest, join them together and take the difference
# lca_remind_reference <- lca_remind %>% filter(scenario == 'Base')
# #lca_remind_reference <- lca_remind %>% filter(scenario == 'KSP90aE')
# 
# #rename value to value_reference
# names(lca_remind_reference)[names(lca_remind_reference) == 'value'] <- 'value_reference'
# 
# #join lca_remind_reference with the lca_remind
# lca_remind <- dplyr::left_join(lca_remind,lca_remind_reference, by=c('region','period','category','end_category','tech','X','RegionCode') )
# 
# lca_remind$value <- lca_remind$value - lca_remind$value_reference
# 
# #rename some of the joining stuff
# names(lca_remind)[names(lca_remind) == 'scenario.x'] <- 'scenario'
# lca_remind$scenario.y <- NULL

#clear names of categories and end_categories
lca_remind$category     <- plyr::mapvalues(lca_remind$category, from = category,
                                       to = category_clear)

lca_remind$end_category <- plyr::mapvalues(lca_remind$end_category, from = end_category,
                                       to = end_category_clear)

#SE$category     <- plyr::mapvalues(SE$category, from = category,
 #                                          to = category_clear)

#SE$end_category <- plyr::mapvalues(SE$end_category, from = end_category,
#                                           to = end_category_clear)
#SE<-SE %>% filter(category %in% category_clear_impact)
lca_remind<-lca_remind %>% filter(category %in% category_clear_impact)

#lca_remind <- lca_remind%>% filter(region=='DEU')
#lca_remind <- lca_remind%>% filter(tech=='SE|Electricity')
#lca_remind <- lca_remind%>% filter(period %in% c('2015','2030','2050'))

#replace NAs with zeros

lca_remind$value[lca_remind$value %in% NA] <- 0
#lca_remind$value_reference[lca_remind$value_reference %in% NA] <- 0
#SE$value[SE$value %in% NA] <- 0



#add some different totals
lca_remind_end_category_tech_total <- lca_remind %>% 
                                      group_by(scenario, RegionCode, period,tech, end_category) %>% 
                                      summarise(value = sum(value))

lca_remind_end_category_total <- lca_remind %>% 
  group_by(scenario, RegionCode, period, end_category) %>% 
  summarise(value = sum(value))

lca_remind_end_category_category_total <- lca_remind %>% 
  group_by(scenario, RegionCode, period, category, end_category) %>% 
  summarise(value = sum(value))

lca_remind_end_category_category_total_region <- lca_remind %>% 
  group_by(scenario, period, category, end_category) %>% 
  summarise(value = sum(value))

lca_iso_end_category_total <- lca_remind %>% 
  group_by(scenario, region, period, end_category) %>% 
  summarise(value = sum(value))

write.csv(lca_iso_end_category_total,paste0(dir_outtab,'lca_indicators_cost.csv'))
#loop over end categories since they have different scales
#loop over sectors?

tech = unique(lca_remind$tech)
end_category = unique(lca_remind$end_category)
category = unique(lca_remind$category)
# add a sector total
#aggregate( df[,11:200], df[,1:10], FUN = sum )

#write results
#write.csv(lca_remind,paste0(dir_outtab,'environmental_indicators_cost.csv'))

# calculate the tech totals

lca_remind_cost_per_tec                          <- lca_remind %>% 
                                                    dplyr::group_by(scenario,tech, RegionCode, period) %>% 
                                                    dplyr::summarise(value = sum(value))



# lca_remind_cost_per_tec_ref  <- lca_remind %>% 
#   dplyr::group_by(scenario,tech, RegionCode, period) %>% 
#   dplyr::summarise(value_reference = sum(value_reference))
# why xxxx
#SE$SE <- SE$SE /54

lca_remind_cost_per_tec <- dplyr::left_join(lca_remind_cost_per_tec,lca_remind_cost_per_tec_ref)


#xxx use here the original SE
SE_gdx$X <- NULL
colnames(SE_gdx) <- c('scenario', 'period', 'RegionCode','tech','SE') 

#### check this, why are some tax 0?
#because they are not available in Base


#replace zeros and NAs with low number
SE_gdx$SE[SE_gdx$SE == 0] <- 10^-9

#from TWyears to EJ to GJ
SE_gdx$SE <- SE_gdx$SE * 31.56 * 10^9


lca_remind_cost_per_tec <- dplyr::left_join(lca_remind_cost_per_tec, SE_gdx)

#10^-9 because we need the tax for GJ
#to GJ
#lca_remind_cost_per_tec$tax <- ((lca_remind_cost_per_tec$value_reference+lca_remind_cost_per_tec$value))/lca_remind_cost_per_tec$SE

lca_remind_cost_per_tec$tax <- ((lca_remind_cost_per_tec$value))/lca_remind_cost_per_tec$SE

out <- lca_remind_cost_per_tec %>% select(c('period', 'RegionCode', 'tech', 'tax'))


out$tax[is.na(out$tax)] <- 0


# that should be redundant?
#take the 2015 values from Base for all other scenarios
out_2015 <- out %>% filter(scenario == 'Base') %>% filter(period == 2015)

out  <- out %>% filter(period != 2015)

for(scen in unique(out$scenario)){
  
  out_2015$scenario <- as.character(scen)
  out <- rbind(out,out_2015)
}


# map here the scenarios to rcps
# improve to take it out of the gdx

out$scenario <- as.character(out$scenario) 

out$scenario[out$scenario == 'Base'] <- 'none'
out$scenario[out$scenario == 'NDC'] <- 'rcp45'
out$scenario[out$scenario == 'NPi'] <- 'rcp45'
out$scenario[out$scenario == 'PkBudg900'] <- 'rcp20'
out$scenario[out$scenario == 'PkBudg1100'] <- 'rcp26'
out$scenario[out$scenario == 'PkBudg1300'] <- 'rcp26'

names(out) <- c("scenario","period","RegionCode","tech","value" )



#interpolate missing values
out <- interpolate_missing_periods(out,method = 'linear' , expand.values=T, period=unique(SE_gdx$period))

# take the world average if under 30% of world average which would probalby mean missing data
out_glo <- out %>% 
  dplyr::group_by(scenario, tech,period) %>% 
  dplyr::summarize(value = mean(value))
out_glo$RegionCode <- 'GLO'
out <- left_join(out, out_glo, by = c("scenario",  "tech", "period"))

out$value <- if_else(out$value.x > out$value.y*0.3, out$value.x,out$value.y)
out <- out %>% select(c("scenario",'period',"RegionCode.x","tech","value" ))
names(out) <- c("scenario",'period',"RegionCode","tech","value" )


# take the 2050 values for everything after because this is missing currently from the lca

out_lt_2050 <- out%>% filter(period <2050)
out_2050 <- out %>% filter(period == 2050)

out_gt_2050 <- out %>% filter(period>2050)

out_gt_2050 <- left_join(out_gt_2050,out_2050,by = c("scenario",  "tech",'RegionCode'))
out_gt_2050 <- out_gt_2050 %>% select(c("scenario",'period.x',"RegionCode","tech","value.y" ))
names(out_gt_2050) <- c("scenario",'period',"RegionCode","tech","value" )

out <- rbind(out_lt_2050,out_2050,out_gt_2050)
  

#this doesnt delete duplicates of only a-d, exclude value colmn for duplicates removal
out <-out[!duplicated(out, by=c("scenario", "period", "RegionCode", "tech")),]


#separatly write out FE and SE
#filter here FE and SE


out_SE <- out %>% filter(tech %in% SE_tech)
out_FE <- out %>% filter(tech %in% FE_tech)

#FE has three techs, seperated here

out_FE <- separate(out_FE,tech,
         c('tech1','tech2','tech3'),
         sep = "_",
         remove = TRUE)

write.table(out_SE,paste0(dir_outtab,'f21_tau_pe2se_tax.cs4r'),row.names = F, col.names = F, quote=FALSE ,sep=',')
write.table(out_FE,paste0(dir_outtab,'f21_tau_fe_tax.cs4r'),row.names = F, col.names = F, quote=FALSE ,sep=',')







if (do.plots) {
  

x<-1
#sectoral category
for(cat in category){
  lca_plot <- lca_remind %>% 
    filter(category == cat)
  
  p <- ggplot(data=lca_plot, aes(x=period, y=value, group = scenario, fill = scenario)) +
    labs(title = paste(cat))+
    geom_bar(stat="identity", position = "dodge")  +
    facet_wrap(~tech) +
    theme(legend.position = 'bottom',axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(filename=paste0(dir_plot_out,"lca_cost_sectoral_category_",cat,".pdf"),p,scale=1,width=25,height=25,unit="cm")
  x<-x+1
}



x<-1

# sector midpoints
for(sector in tech){
lca_plot <- lca_remind %>%
  #filter(period == as.integer('2050')) %>%
  filter(tech == sector)

p <- ggplot(data=lca_plot, aes(x=period, y=value, group = scenario, fill = scenario)) +
     labs(title = paste(sector))+
     geom_bar(stat="identity", position = "dodge")  +
     facet_wrap(~category) +
     theme(legend.position = 'bottom',axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename=paste0(dir_plot_out,"lca_cost_sectoral_",x,".pdf"),p,scale=1,width=25,height=25,unit="cm")
x<-x+1
}

x<-1

#sector endpoints
for(sector in tech){
  lca_plot <- lca_remind_end_category_tech_total %>% 
    #filter(period == as.integer('2050')) %>%
    filter(tech == sector)
  
  p <- ggplot(data=lca_plot, aes(x=period, y=value, group = scenario, fill = scenario)) +
    labs(title = paste(sector))+
    geom_bar(stat="identity", position = "dodge")  +
    facet_wrap(~end_category) +
    theme(legend.position = 'bottom',axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(filename=paste0(dir_plot_out,"lca_cost_end_category_total_",x,".pdf"),p,scale=1,width=25,height=25,unit="cm")
  x<-x+1
}
x<-1
#total category
for(end in end_category){
  lca_plot <- lca_remind_end_category_category_total %>% 
    filter(end_category == end)
  
  p <- ggplot(data=lca_plot, aes(x=period, y=value, group = scenario, fill = scenario)) +
    labs(title = paste(end))+
    geom_bar(stat="identity", position = "dodge")  +
    facet_wrap(~category) +
    theme(legend.position = 'bottom',axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(filename=paste0(dir_plot_out,"lca_cost_end_category_category_",end,".pdf"),p,scale=1,width=25,height=25,unit="cm")
  x<-x+1
}

#total end_category

  lca_plot <- lca_remind_end_category_total %>% filter(period < 2100)
  
 p<- ggplot(data=lca_plot, aes(x=period, y=value, group = scenario, fill = scenario)) +
    geom_bar(stat="identity", position = "dodge")  +
    facet_wrap(~end_category) +
    theme(legend.position = 'bottom',axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(filename=paste0(dir_plot_out,"lca_cost_remind_end_category_total.pdf"),p,scale=1,width=25,height=25,unit="cm")
  
#filter end_category policy cost for one policy cost  
  lca_plot <- lca_remind_end_category_category_total_region%>% filter(period < 2100)
  lca_plot <- lca_plot %>% filter(end_category %in% c(end_category_clear[1:3],end_category_clear[7:7]))%>% filter(end_category !=  'Ressource Depletion [US$2005]')
  lca_plot$value <- lca_plot$value*-1
  lca_plot_total <- lca_plot %>% 
    group_by(scenario, period) %>% 
    summarise(total = sum(value))
  
  
p<-ggplot(data=lca_plot, aes(x=period, y=value)) +
    geom_col(position = "stack",aes(group = category, fill = as.factor(category)))  +
    geom_point(data = lca_plot_total,aes(x=period, y=total))+
    facet_grid(~scenario) 
    theme(legend.position = 'bottom',axis.text.x = element_text(angle = 90, hjust = 1))
   ggsave(filename=paste0(dir_plot_out,"lca_cost_remind_end_category_total_categorie.pdf"),p,scale=1,width=25,height=25,unit="cm")

   
   
   #filter end_category policy cost for one policy cost, region wrap
   lca_plot <- lca_remind_end_category_category_total%>% filter(period < 2100)
   lca_plot <- lca_plot %>% filter(end_category %in% c(end_category_clear[1:3],end_category_clear[7:7]))%>% filter(end_category !=  'Ressource Depletion [US$2005]')
   lca_plot$value <- lca_plot$value*-1
   lca_plot_total <- lca_plot %>% 
     group_by(scenario, period,RegionCode) %>% 
     summarise(total = sum(value))
   
 p <- ggplot(data=lca_plot, aes(x=period, y=value)) +
     geom_col(position = "stack",aes(group = category, fill = as.factor(category)))  +
     geom_point(data = lca_plot_total,aes(x=period, y=total))+
     facet_wrap(~scenario~RegionCode) +
     theme(legend.position = 'bottom',axis.text.x = element_text(angle = 90, hjust = 1))
   ggsave(filename=paste0(dir_plot_out,"lca_cost_remind_end_category_total_categorie_region.pdf"),p,scale=1,width=40,height=40,unit="cm")
}

