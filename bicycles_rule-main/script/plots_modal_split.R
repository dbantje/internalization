library(ggplot2)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(cowplot)

#make nice graphs, windrose?

df <- read.csv("data/Budg1100_ElecEra_2020-09-17_17.17.49/REMIND_generic_Budg1100_ElecEra_withoutPlus.mif",sep=";")

c1 <- c("FE","FE|Transport|Freight","FE|Transport|Pass")
c2 <- c("ES|Transport|Pass","ES|Transport|Pass|Rail","ES|Transport|Pass|Aviation|Domestic","ES|Transport|Pass|Aviation","ES|Transport|Pass|Road")
c3 <- c("ES|Transport|Pass|Road","ES|Transport|Pass|Road|LDV","ES|Transport|Pass|Road|LDV|BEV","ES|Transport|Pass|Road|LDV|FCEV","ES|Transport|Pass|Road|LDV|Gases","ES|Transport|Pass|Road|LDV|Hybrid Electric","ES|Transport|Pass|Road|LDV|Hybrid Liquids","ES|Transport|Pass|Road|LDV|Liquids","ES|Transport|Pass|Road|Non-Motorized")



df <- df %>% filter(Region=="World")

df1 <- df %>% filter(Variable %in%  c1)
df1$group <- "Final energy"
df2 <- df %>% filter(Variable %in%  c2)
df2$group <- "Passenger Transport"
df3 <- df %>% filter(Variable %in%  c3)
df3$group <- "Modal split"

df <- rbind(df1,df2,df3)

df <- melt(df, c(id.vars=colnames(df)[1:5],"group"))


df$value <- as.numeric(df$value)



scenarios <- unique(df$Scenario)

#stacked FE

df <- df %>% filter(variable %in%  c("X2020","X2030"))
p_FE <- ggplot(df,aes(y=value,x=Variable, fill = factor(Variable), group=Variable)) +
        geom_bar(data = df %>% filter(Variable %in%  c1[1])%>% filter(Scenario %in% scenarios[1]),
                           aes(x=Variable, y=value, fill = factor(Variable), group=Scenario),position = "dodge2",stat = 'identity') +
        geom_bar(data = df  %>% filter(Variable %in%  c1[2:3])%>% filter(Scenario %in% scenarios[1]),
                 aes(x=Scenario, y=value, fill = factor(Variable)),position = "stack",stat = 'identity') +
        ylab("EJ/yr")+
        theme(legend.position = 'bottom',axis.text.x = element_text(angle = 30, hjust = 1)) +
        facet_wrap(~variable)

ggsave("plot/transport_FE.png",p_FE, width = 12, height = 12)
ggsave("plot/transport_FE.pdf",p_FE, width = 12, height = 12)
df <- df %>% filter(variable %in%  c("X2030"))
df_bike <- df
df_bike$Scenario <- paste(unique(df$Scenario),"_urban_LDV_80%")

bike_pkm <-  df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV","value"] * 0.47263520008

df_bike[df_bike$Variable == "ES|Transport|Pass|Road","value"]                      <- df_bike[df_bike$Variable == "ES|Transport|Pass|Road","value"]   -bike_pkm
df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV","value"]                  <- df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV","value"] * (1-0.47263520008)
df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV|BEV","value"]              <- df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV|BEV","value"] * (1-0.47263520008)
df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV|FCEV","value"]             <- df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV|FCEV","value"] * (1-0.47263520008)
df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV|Gases","value"]            <- df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV|Gases","value"] * (1-0.47263520008)
df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV|Hybrid Electric","value"]  <- df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV|Hybrid Electric","value"] * (1-0.47263520008)
df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV|Hybrid Liquids","value"]   <- df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV|Hybrid Liquids","value"] * (1-0.47263520008)
df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV|Liquids","value"]          <- df_bike[df_bike$Variable == "ES|Transport|Pass|Road|LDV|Liquids","value"] * (1-0.47263520008)
df_bike[df_bike$Variable == "ES|Transport|Pass|Road|Non-Motorized","value"]        <- df_bike[df_bike$Variable == "ES|Transport|Pass|Road|Non-Motorized","value"] +  bike_pkm

df_bike[df_bike$Variable == "ES|Transport|Pass|Road|Non-Motorized","Variable"]     <- "ES|Transport|Pass|Active Transport"
df[df$Variable == "ES|Transport|Pass|Road|Non-Motorized","Variable"]     <- "ES|Transport|Pass|Active Transport"

df <- rbind(df,df_bike)



p_transport <- ggplot(df,aes(y=value,x=Variable, fill = factor(Variable), group=Scenario)) +
   geom_bar(data = df  %>% filter(Variable %in%  c2[2:5]) ,
           aes(x=Scenario, y=value, fill = factor(Variable)),position = "stack",stat = 'identity') +
  geom_bar(data = df  %>% filter(Variable %in%  c(c3[2:9], "ES|Transport|Pass|Active Transport")) ,
           aes(x=Scenario, y=value, fill = factor(Variable)),position = "dodge2",stat = 'identity') +
  ylab("bn pkm/yr")+
  theme(legend.position = 'bottom',axis.text.x = element_text(angle = 30, hjust = 1)) 


ggsave("plot/transport_mode.png",p_transport, width = 12, height = 12)
ggsave("plot/transport_mode.pdf",p_transport, width = 12, height = 12)

p_transport_total <- plot_grid(p_FE,p_transport , align = "v",axis = "b", nrow = 1,ncol=2)

ggsave("plot/transport_mode_total.png",p_transport_total, width = 12, height = 12)
ggsave("plot/transport_mode_total.pdf",p_transport_total, width = 12, height = 12)

write.csv(df,"transport_mode_total_data.csv")
