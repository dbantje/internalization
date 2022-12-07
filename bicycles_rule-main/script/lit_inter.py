#scaling factors
pop_scaling_2010_EU = remind.query('Variable == ("Population") and Region.str.contains("World")', engine='python')[str(year)]/remind.query('Variable == ("Population") and Region.str.contains("EUR")', engine='python')['2010'].iloc[0]

pop_scaling_2020 = remind.query('Variable == ("Population") and Region.str.contains("World")', engine='python')[str(year)]/remind.query('Variable == ("Population") and Region.str.contains("World")', engine='python')['2020']

pop_scaling_2015 = remind.query('Variable == ("Population") and Region.str.contains("World")', engine='python')[str(year)]/remind.query('Variable == ("Population") and Region.str.contains("World")', engine='python')['2015']
pop_scaling_2040 = remind.query('Variable == ("Population") and Region.str.contains("World")', engine='python')[str(year)]/remind.query('Variable == ("Population") and Region.str.contains("World")', engine='python')['2040']

#change this here to 
vkm_scaling_2020 = remind.query('Variable == ("ES|Transport|VKM|Pass|Road|LDV") and Region.str.contains("World")', engine='python')[str(year)]/remind.query('Variable == ("ES|Transport|VKM|Pass|Road|LDV") and Region.str.contains("World")', engine='python')['2020']

ownership_cost_scaling = remind.query('Variable == ("GDP|per capita|MER") and Region.str.contains("World")', engine='python')[str(year)]/remind.query('Variable == ("GDP|per capita|MER") and Region.str.contains("EUR")', engine='python')['2020'].iloc[0]
ownership_cost = remind.query('Variable == ("Price|Energy Service|Transport LDV") and Region.str.contains("World")', engine='python')[str(year)] * remind.query('Variable == ("FE|Transport|Pass|Road|LDV") and Region.str.contains("World")', engine='python')[str(year)].iloc[0] * 1000000000
congestion = remind.query('Variable == ("GDP|MER") and Region.str.contains("World")', engine='python')[str(year)] 

tire_microplastic_scaling = ( (remind.query('Variable == ("ES|Transport|VKM|Pass|Road") and Region.str.contains("World")', engine='python')[str(year)] + remind.query('Variable == ("ES|Transport|VKM|Freight|Road") and Region.str.contains("World")', engine='python')[str(year)].iloc[0]) / (remind.query('Variable == ("ES|Transport|VKM|Pass|Road") and Region.str.contains("World")', engine='python')['2015'] + remind.query('Variable == ("ES|Transport|VKM|Freight|Road") and Region.str.contains("World")', engine='python')['2015'].iloc[0]))

tire_microplastic_scaling_1 = remind.query('Variable == ("ES|Transport|VKM|Pass|Road|LDV") and Region.str.contains("World")', engine='python')[str(year)].iloc[0]/ (remind.query('Variable == ("ES|Transport|VKM|Pass|Road") and Region.str.contains("World")', engine='python')[str(year)] + remind.query('Variable == ("ES|Transport|VKM|Freight|Road") and Region.str.contains("World")', engine='python')[str(year)].iloc[0])

Physical_activity_health_impact_GBD_DALY = literature.query('variable == ("Physical_activity_health_impact_GBD") and unit.str.contains("DALY")', engine='python')
Physical_activity_health_impact_GBD_DALY['value'] =  literature.query('variable == ("Physical_activity_health_impact_GBD") and unit.str.contains("DALY")', engine='python')['value_lit'][0] * pop_scaling_2015.iloc[0,]

Physical_activity_health_impact_GBD_premature_death = literature.query('variable == ("Physical_activity_health_impact_GBD") and unit.str.contains("premature death")', engine='python')
Physical_activity_health_impact_GBD_premature_death['value'] =  literature.query('variable == ("Physical_activity_health_impact_GBD") and unit.str.contains("premature death")', engine='python')['value_lit'] * pop_scaling_2015.iloc[0,]

Physical_activity_health_care_cost_WHO = literature.query('variable == ("Physical_activity_health_care_cost_WHO") and unit.str.contains("$")', engine='python')
Physical_activity_health_care_cost_WHO['value'] =  literature.query('variable == ("Physical_activity_health_care_cost_WHO") and unit.str.contains("$")', engine='python')['value_lit'] * pop_scaling_2015.iloc[0,]

Physical_activity_health_impact_lancet = literature.query('variable == ("Physical_activity_health_impact_lancet") and unit.str.contains("premature death")', engine='python')
Physical_activity_health_impact_lancet['value'] =  literature.query('variable == ("Physical_activity_health_impact_lancet") and unit.str.contains("premature death")', engine='python')['value_lit'] * pop_scaling_2040.iloc[0,]

Road_accidents_DALY = literature.query('variable == ("Road_accidents") and unit.str.contains("DALY")', engine='python')
Road_accidents_DALY['value'] =  literature.query('variable == ("Road_accidents") and unit.str.contains("DALY")', engine='python')['value_lit'] * pop_scaling_2020.iloc[0,] * vkm_scaling_2020.iloc[0,]

Road_accidents_premature_deaths = literature.query('variable == ("Road_accidents") and unit.str.contains("premature death")', engine='python')
Road_accidents_premature_deaths['value'] =  literature.query('variable == ("Road_accidents") and unit.str.contains("premature death")', engine='python')['value_lit'] * pop_scaling_2020.iloc[0,]* vkm_scaling_2020.iloc[0,]

Road_noise = literature.query('variable == ("Road_noise") and unit.str.contains("DALY")', engine='python')
Road_noise['value'] =  literature.query('variable == ("Road_noise") and unit.str.contains("DALY")', engine='python')['value_lit'] * pop_scaling_2010_EU.iloc[0,]* vkm_scaling_2020.iloc[0,]

Understimated_additional_owner_cost = literature.query('variable == ("Understimated_additional_owner_cost") and unit.str.contains("$")', engine='python')
Understimated_additional_owner_cost['value'] =  literature.query('variable == ("Understimated_additional_owner_cost") and unit.str.contains("$")', engine='python')['value_lit'] * ownership_cost_scaling.iloc[0,] * ownership_cost.iloc[0,]

Congestion_fuel_time_waste = literature.query('variable == ("Congestion_fuel_time_waste") and unit.str.contains("$")', engine='python')
Congestion_fuel_time_waste['value'] =  literature.query('variable == ("Congestion_fuel_time_waste") and unit.str.contains("$")', engine='python')['value_lit'] * congestion.iloc[0,] * 1000000000

tire_microplastic = literature.query('variable == ("Platic_emission_tire") and unit.str.contains("t")', engine='python')

tire_microplastic['value'] = tire_microplastic_scaling.iloc[0] * tire_microplastic_scaling_1.iloc[0]  *  literature.query('variable == ("Platic_emission_tire") and unit.str.contains("t")', engine='python')['value_lit'].iloc[0]

literature = pd.concat([Physical_activity_health_impact_GBD_DALY,
                                   Physical_activity_health_impact_GBD_premature_death,
                                   Physical_activity_health_care_cost_WHO,
                                   Physical_activity_health_impact_lancet,
                                   Road_accidents_DALY,
                                   Road_accidents_premature_deaths,
                                   Road_noise,
                                   Understimated_additional_owner_cost,
                                   Congestion_fuel_time_waste,
								   tire_microplastic])  
literature['scenario'] = "internalization_modes_total"
literature = literature[['scenario', 'variable', 'unit','value']]


literature['variable'] = literature['variable'] + '_' + literature['unit']
literature = literature[['scenario', 'variable','value']]

literature=literature.reset_index(drop=True)
literature=literature.pivot(index='variable', columns=['scenario'], values='value')
