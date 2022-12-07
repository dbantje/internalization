# noise, extrapolate from EU-33 https://www.eea.europa.eu/publications/health-risks-caused-by-environmental
# DALYS 2017 (middle of ranges)

#Road 	High annoyance 275500
#High sleep disturbance 194500
#Ischaemic heart disease 131000


#Rail 	High annoyance 52500
#High sleep disturbance 77000
#Ischaemic heart disease 21100


#Aircraft 	High annoyance 16500
#High sleep disturbance 11000
#Ischaemic heart disease 2425

edge_internalization_mapping_noise           = {"Compact Car": "Road",
                                               'Large Car': "Road",
                                               'Large Car and SUV': "Road",
                                               'Light Truck and SUV': "Road",
                                               'Midsize Car': "Road",
                                               'Mini Car': "Road",
                                               'Moped': 'Road',
                                               'Motorcycle (50-250cc)': 'Road',
                                               'Motorcycle (>250cc)': 'Road',
                                               'Subcompact Car': "Road",
                                               'Van': "Road",
                                               'Bus_tmp_vehicletype': "Road",
                                               'Freight Rail_tmp_vehicletype': "Rail",
                                               'HSR_tmp_vehicletype': "Rail",
                                               'Passenger Rail_tmp_vehicletype': "Rail",
                                               'Truck (0-3.5t)': "Road",
                                               'Truck (18t)': "Road",
                                               'Truck (26t)': "Road",
                                               'Truck (40t)': "Road",
                                               'Truck (7.5t)': "Road",
                                               'Domestic Aviation_tmp_vehicletype': "Aircraft" ,
                                               'International Aviation_tmp_vehicletype': "Aircraft"}

df_noise =  pd.read_csv(noise_eea_paths, sep = ';'   , encoding='latin-1')

REMIND = ["ES|Transport|VKM|Pass|Aviation|Domestic",
          "ES|Transport|VKM|Pass|Road",
          "ES|Transport|VKM|Pass|Rail"]

#pedestrian is assigned to cars, its therefore missing in the mapping
noise_REMIND_mapping                        = {"ES|Transport|VKM|Pass|Aviation|Domestic": "Aircraft",
                                               "ES|Transport|VKM|Pass|Road": "Road",
                                               "ES|Transport|VKM|Pass|Rail": "Rail"}

#filter the remind results for the activities abvove
remind_noise      = remind.loc[(remind['Variable'].isin(REMIND))]

#map
remind_noise['Variable'] = remind_noise['Variable'].map(noise_REMIND_mapping)
remind_noise = pd.melt(remind_noise, id_vars=['Model','Scenario','Region' ,'Variable','Unit'], var_name='Period', value_name='value',ignore_index=True)

#delete 
remind_noise    = remind_noise.loc[(~remind_noise['Period'].isin(['Unnamed: 24']))]
remind_noise['Period'] = remind_noise['Period'].astype(int)

#merge
df_noise = df_noise.merge(remind_noise,  how='left',on = ['Region','Period','Variable'])
df_noise['value'] = df_noise['Value'] / (df_noise['value'])
df_noise['Unit']  = 'Noise DALY/bn vkm'
# scale all years and regions with this noise result for EUR

#join this with remind results regardless of region and period
df_noise = df_noise.drop(['Unit_x','Value','Unit_y'],axis=1)
df_noise = df_noise.merge(remind_noise,  how='right',on = ['Variable'])
df_noise["value"] = df_noise["value_y"]

df_noise =  df_noise[['Region_y','Period_y','Scenario_y','Variable','Unit_x','value']]
df_noise = df_noise.rename(columns={"Region_y": "Region","Period_y":"Period","Scenario_y": "Scenario","Unit_x":"Unit"})

