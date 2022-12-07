library(ggplot2)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(tidyverse)

#make nice graphs, windrose?

df <- read.csv("results/data_1.csv",sep=";")

#add a monetized health
health_monetized <- df %>% filter(end_point_unit == 'DALY') 
health_monetized$end <- health_monetized$end * 72466.4

health_monetized$end_point_unit <- '$'


#add a monetized of environment
env_monetized <- df %>% filter(end_point_unit == 'species.yr') 
env_monetized$end <- env_monetized$end * 11148648.72

env_monetized$end_point_unit <- '$'

df <- df %>% filter(end_point_unit != 'species.yr') 


df <- rbind(df,health_monetized,env_monetized)

#plot
p<- list()
for (x in unique(df$end_point)) {
  
  df_tmp <-  df %>% filter(end_point == x)
  p[[x]]<- ggplot(df_tmp, aes(y=end,x=indicator, fill = factor(scenario))) +
    geom_bar(position = "dodge2",stat = 'identity') +
    theme(legend.position = 'bottom',axis.text.x = element_text(angle = 30, hjust = 1)) +
    facet_wrap(~ end_point + end_point_unit,scales = 'free')
  
}
p <- do.call(grid.arrange,p)


p_climate_total <-p[[1]]$climate


ggsave("plot/impact_total_bike_share_0.47263520008_2030_lit.png",p, width = 13, height = 20)
ggsave("plot/impact_total_bike_share_0.47263520008_2030_lit.pdf",p, width = 13, height = 20)




## per pkm lit impacts

#read per tec lca results
pkm <- read.csv("results/impact_pkm_2030.csv",sep=",")
#read remind vkm ldv
remind <- read.csv("data/Budg1100_ElecEra_2020-09-17_17.17.49/REMIND_generic_Budg1100_ElecEra.mif",sep=";")
remind <- remind %>% filter(Region == 'World')

ldv_road <- remind %>% filter(Variable == "ES|Transport|Pass|Road|LDV")

ldv_road <- as.numeric(ldv_road$X2030)

#distribute lit impacts on the ldv
df$end <- df$end/(as.numeric(ldv_road)*1000000000)
df <- df%>% filter(scenario == 'ldv_total')



#wide to long
pkm <- melt(pkm, id.vars=c("indicator"))
pkm$"mode" <- pkm$variable
pkm$variable <- NULL

#replace the original pkm and calc the end
pkm$pkm <-pkm$value
pkm$value <-NULL

df$gasoline_car <- df$end
df$battery_electric_car <- df$end
df$bike <- df$end*0
df$ebike <- df$end*0

df$end <- NULL

#wide to long
df <- melt(df, id.vars=c("scenario","indicator","category","unit","mid","end_point_factor","end_point","end_point_unit"))

df$mode <- df$variable
df$variable <- NULL


air_poll_share <- pkm %>% filter(indicator == 'particulate matter formation')

pkm <- left_join(pkm,df, by = c('indicator','mode'))
pkm$value <- pkm$pkm * pkm$end_point_factor
pkm <- drop_na(pkm)

pkm <- pkm %>% select(scenario, indicator,end_point, end_point_unit,mode, value)
df <-  df  %>% select(scenario, indicator,end_point, end_point_unit,mode, value)

indicators <- unique(pkm$indicator)


'%!in%' <- function(x,y)!('%in%'(x,y))
df <- df %>% filter(indicator %!in% indicators)

#take here the share of particulate matter between the BEV and gasoline car modes as proxy for air pollution related impacts


#take average particulate matter formation from all modes and calc share per mode
#take all air containing indicators and apply this share
air_poll_share$mean <- mean(air_poll_share$pkm)
air_poll_share$pkm  <- air_poll_share$pkm/air_poll_share$mean

air_poll_share <- air_poll_share %>% select(mode,pkm)
#########
air_poll_df  <- df %>% filter(str_detect(indicator, "air"))


air_poll_df <- left_join(air_poll_df, air_poll_share, by="mode")
air_poll_df$value <- air_poll_df$value * air_poll_df$pkm
air_poll_df$pkm <- NULL

df  <- df %>% filter(!str_detect(indicator, "air"))

df  <- rbind(df,air_poll_df )

# 
# #monetize env
pkm$value[which(pkm$end_point == 'environment')] <- pkm$value[which(pkm$end_point == 'environment')] * 11148648.72

df <- rbind(pkm,df)

#filter out air pollution from economic
df$end_point[which(df$indicator %in% "air_pollution_cost")] <- 'health'
df$end_point[which(df$indicator %in% "air_pollution_direct_cost")] <- 'health'
df$end_point[which(df$indicator %in% "air_pollution_indirect_cost")] <- 'health'

df$end_point[which(df$indicator %in% "Physical_activity_health_care_cost_WHO")] <- 'health'

df <-df[which( df$end_point_unit %!in% 'DALY'),]
df <-df[which(df$end_point_unit %!in% 'premature death'),]


df$mode <- as.character(df$mode)


df <- df %>% distinct(indicator, end_point, end_point_unit,mode, .keep_all = TRUE)

#colour BEV and gasoline car with similar colours but different hues
cols<-c("battery_electric_car" = "#E69F00","bike"= "#1C9371","ebike"= "#a3d71c","gasoline_car"= "#e35d29","car"= "#e97f56")

# limit to climate, resources, maybe tire plastic
df[df$indicator == "Platic_emission_tire","end_point"]                      <- "ressources"

    df_tmp <-  df %>% filter(end_point == 'climate')
    p_climate<- ggplot(df_tmp, aes(y=value,x=indicator, fill = factor(mode))) +
               geom_bar(position = "dodge2",stat = 'identity') +
               xlab(NULL)+
               ylab(NULL)+
               theme_minimal()+ 
               theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
               scale_fill_manual(values=cols)+
               theme(legend.position = "bottom")+
               theme(legend.title=element_blank())+
               labs(title="Climate \n [kg CO2 eq]") 
    
    df_tmp <-  df %>% filter(end_point == 'ressources')
    p_ressources<- ggplot(df_tmp, aes(y=value,x=indicator, fill = factor(mode))) +
      geom_bar(position = "dodge2",stat = 'identity') +
      xlab(NULL)+
      ylab(NULL)+
      theme_minimal()+ 
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      scale_fill_manual(values=cols)+
      theme(legend.position = "bottom") +
      facet_wrap(~ end_point + end_point_unit,scales = 'free',ncol=4)
    df_tmp <-  df %>% filter(indicator == 'Platic_emission_tire')
    # p_tire<- ggplot(df_tmp, aes(y=value,x=indicator, fill = factor(mode))) +
    #   geom_bar(position = "dodge2",stat = 'identity') +
    #   xlab(NULL)+
    #   ylab(NULL)+
    #   theme_minimal()+ 
    #   theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    #   scale_fill_manual(values=cols)+
    #   theme(legend.position = "none") +
    #   theme(legend.title=element_blank())
  




# do rest end_points
# outlines like colors above
# seperate legend in end categories


end <- c("economic",    "environment", "health")
  
  df_tmp <-  df %>% filter(end_point %in% end) %>% filter(indicator %!in% 'air_pollution') %>% filter(indicator %!in% 'air_pollution_cost') 
  
 # add a end_point total?
  # total <- df_tmp %>% 
  #   group_by(end_point,mode) %>% 
  #   summarise(value = sum(value))
  # total$"scenario" <- unique(df_tmp$scenario)
  # total$"end_point_unit" <- unique(df_tmp$end_point_unit)
  # total$"indicator" <- total$end_point
  # total$end_point <- 'total'
  # 
  # df_tmp <- rbind(df_tmp,total)

  
  p_cost <- ggplot(df_tmp, aes(y=value,x=mode, fill = factor(indicator), group = mode)) +
    #geom_bar( aes(y=value,x=mode, fill = factor(end_point), group = end_point), position = "stack",stat = 'identity',alpha =1, width = 1.5) +
    geom_bar(position = "stack",stat = 'identity') +
    
    
    xlab(NULL)+
    ylab(NULL)+
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    theme(legend.position = "bottom") +
    theme(legend.title=element_blank())+
    #scale_fill_manual(values=cols)+
    facet_wrap(~ end_point ,scales = 'free',ncol = 3)
  
  # p_cost_total <- ggplot(df_tmp, aes(y=value,x=mode, fill = factor(indicator), group = mode)) +
  #   geom_bar( aes(y=value,x=mode, fill = factor(end_point), group = end_point), position = "stack",stat = 'identity',alpha =1, width = 1.5) +
  #   
  #  
  #   
  #   xlab(NULL)+
  #   ylab(NULL)+
  #   theme_minimal()+ 
  #   theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  #   theme(legend.position = "bottom")+
  #   theme(legend.title=element_blank())+
  #   labs(title="Cost \n [$]")
  #   #scale_fill_manual(values=cols)+



p  <- plot_grid(p_cost,p_ressources , align = "v", nrow = 2, rel_heights = c(0.8,0.5))

#plot



ggsave("plot/impact_pkm_bike_share_0.47263520008_2030_lit.png",p, width = 10, height =8)
ggsave("plot/impact_pkm_bike_share_0.47263520008_2030_lit.pdf",p, width = 10, height = 8)



################################################################

# save the total CO2 and pkm CO2 emissions here

p  <- plot_grid(p_climate_total, p_climate , align = "v", nrow = 1)

ggsave("plot/climate.png",p, width = 6, height = 6)
ggsave("plot/climate.pdf",p, width = 6, height = 6)

###############################################################

#difference plot
df <- read.csv("results/data_2.csv",sep=";")
df$end <- df$ldv_total - df$ldv_total_bike.share.0.47263520008

#add a monetized health
health_monetized <- df %>% filter(end_point_unit == 'DALY') 
health_monetized$end <- health_monetized$end * 72466.4

health_monetized$end_point_unit <- '$'

df <- rbind(df,health_monetized)

ggplot(df, aes(y=end,x=indicator, fill = 'red', group = end)) +
  geom_bar(position = "dodge2",stat = 'identity') +
  theme(legend.position = 'bottom',axis.text.x = element_text(angle = 30, hjust = 1)) +
  facet_wrap(~ category + end_point_unit,scales = 'free')

ggsave("plot/impact_total_bike_share_0.47263520008_2030_lit_diff.png", width = 12, height = 20)
ggsave("plot/impact_total_bike_share_0.47263520008_2030_lit_diff.pdf", width = 12, height = 20)

#check whats up with ozone depletion
#write
#monetize


