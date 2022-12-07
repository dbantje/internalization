#read air pollution results
air_pollution = pd.read_csv(air_pollution_report, sep=';')
air_pollution = air_pollution[air_pollution["year.1"] == year]

#scale with total road transport
air_pollution_scaling = (1-bike_share/car_load_factor)*  remind.query('Variable == ("ES|Transport|VKM|Pass|Road|LDV") and Region.str.contains("World")', engine='python')[str(year)].iloc[0]/ (remind.query('Variable == ("ES|Transport|VKM|Pass|Road") and Region.str.contains("World")', engine='python')[str(year)] + remind.query('Variable == ("ES|Transport|VKM|Freight|Road") and Region.str.contains("World")', engine='python')[str(year)].iloc[0])

scenario = "ldv_total_bike share "+str(bike_share)


air_pollution_dalys = air_pollution['daly_total'].sum() * air_pollution_scaling.iloc[0]
air_pollution_premature_deaths = air_pollution['data_df_mort_pm.all_mort'].sum()* air_pollution_scaling.iloc[0]
air_pollution_indirect_cost = air_pollution['indirect_daly_cost'].sum()* air_pollution_scaling.iloc[0]
air_pollution_direct_cost = air_pollution['direct_daly_cost'].sum()* air_pollution_scaling.iloc[0]
air_pollution_cost = air_pollution['daly_cost'].sum()* air_pollution_scaling.iloc[0]

d = {'scenario': [scenario, scenario, scenario, scenario, scenario],
     'variable': ['air_pollution_dalys', 'air_pollution_premature_deaths','air_pollution_indirect_cost','air_pollution_direct_cost','air_pollution_cost'],
	 'value': [air_pollution_dalys,air_pollution_premature_deaths,air_pollution_indirect_cost,air_pollution_direct_cost,air_pollution_cost]}

air_pollution_res_bike_share = pd.DataFrame(data=d)

##the same for non bike share


air_pollution_scaling = remind.query('Variable == ("ES|Transport|VKM|Pass|Road|LDV") and Region.str.contains("World")', engine='python')[str(year)].iloc[0]/ (remind.query('Variable == ("ES|Transport|VKM|Pass|Road") and Region.str.contains("World")', engine='python')[str(year)] + remind.query('Variable == ("ES|Transport|VKM|Freight|Road") and Region.str.contains("World")', engine='python')[str(year)].iloc[0])

scenario = "ldv_total"

air_pollution_dalys = air_pollution['daly_total'].sum()* air_pollution_scaling.iloc[0]
air_pollution_premature_deaths = air_pollution['data_df_mort_pm.all_mort'].sum()* air_pollution_scaling.iloc[0]
air_pollution_indirect_cost = air_pollution['indirect_daly_cost'].sum()* air_pollution_scaling.iloc[0]
air_pollution_direct_cost = air_pollution['direct_daly_cost'].sum()* air_pollution_scaling.iloc[0]
air_pollution_cost = air_pollution['daly_cost'].sum()* air_pollution_scaling.iloc[0]

d = {'scenario': [scenario, scenario, scenario, scenario, scenario],
     'variable': ['air_pollution_dalys', 'air_pollution_premature_deaths','air_pollution_indirect_cost','air_pollution_direct_cost','air_pollution_cost'],
	 'value': [air_pollution_dalys,air_pollution_premature_deaths,air_pollution_indirect_cost,air_pollution_direct_cost,air_pollution_cost]}

air_pollution_res = pd.DataFrame(data=d)



air_pollution_res= pd.concat([air_pollution_res,air_pollution_res_bike_share])



air_pollution_res=air_pollution_res.reset_index(drop=True)
air_pollution_res=air_pollution_res.pivot(index='variable', columns=['scenario'], values='value')
