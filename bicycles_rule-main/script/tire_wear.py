# tire wear, just put base emissions on road based transport
# 2.9 Mt in 2014 from road based transport

edge_internalization_mapping_tire            = {'Large Car': "Plastic emissions tires",
                                               'Large Car and SUV': "Plastic emissions tires",
                                               'Light Truck and SUV': "Plastic emissions tires",
                                               'Midsize Car': "Plastic emissions tires",
                                               'Mini Car': "Plastic emissions tires",
                                               'Moped': 'Plastic emissions tires',
                                               'Motorcycle (50-250cc)': 'Plastic emissions tires',
                                               'Motorcycle (>250cc)': 'Plastic emissions tires',
                                               'Subcompact Car': "Plastic emissions tires",
                                               'Van': 'Plastic emissions tires',
                                               'Bus_tmp_vehicletype': 'Plastic emissions tires',
                                               'Truck (0-3.5t)': 'Plastic emissions tires',
                                               'Truck (18t)':    'Plastic emissions tires',
                                               'Truck (26t)':    'Plastic emissions tires',
                                               'Truck (40t)':    'Plastic emissions tires',
                                               'Truck (7.5t)':   'Plastic emissions tires'}

# filter for road based 

# map REMIND to GBD
REMIND = ["ES|Transport|VKM|Pass|Road"]

#filter the remind results for the activities above
remind_tire       = remind.loc[(remind['Variable'].isin(REMIND)& (remind['Region'] == 'World'))]
remind_tire[['2015']] = 2.9/ remind_tire[['2015']] 

#only select 2015 and other columns
remind_tire[['Unit']] = 'Mt/bn vkm'

remind_tire = remind_tire[['Model','Scenario','Region','Variable','Unit','2015']]
remind_tire = remind_tire.rename(columns={"2015":"value"})

