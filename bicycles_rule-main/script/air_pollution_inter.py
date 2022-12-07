#scale the air pollution results to the internalization modes

air_pollution = pd.read_csv(air_pollution_report, sep=';')

air_pollution['year.1']= air_pollution['year.1'].replace(2015,2020)


# air pollution results of the mobility sector
air_pollution = air_pollution[air_pollution["year.1"].isin(years)]

#filter for the relevant scaling remind activities
remind_air_pollution_scaling_variables = ["ES|Transport|VKM|Pass|Road|LDV",
"ES|Transport|VKM|Pass|Road",
"ES|Transport|VKM|Pass|Road|Bus",
"ES|Transport|VKM|Freight|Road"]

air_pollution_variables = ["year.1","iso.iso.y","data_df_mort_pm.all_mort","daly_total","indirect_daly_cost","direct_daly_cost"]

remind_air_pollution = remind[remind['Variable'].isin(remind_air_pollution_scaling_variables)]
remind_air_pollution = remind_air_pollution[['Model','Scenario','Region','Variable','Unit',str(year)]]
remind_air_pollution = remind_air_pollution.pivot(index=['Model','Scenario','Region'], columns='Variable', values=str(year))


#from long to wide with region
#share of considered modes vs the total from the air pollution model
remind_air_pollution['air_pollution_scaling'] = (remind_air_pollution["ES|Transport|VKM|Pass|Road|LDV"] + remind_air_pollution["ES|Transport|VKM|Pass|Road|Bus"]) / ((remind_air_pollution["ES|Transport|VKM|Pass|Road"] + remind_air_pollution["ES|Transport|VKM|Freight|Road"]))

#down to one VKM
remind_air_pollution['air_pollution_scaling'] =  remind_air_pollution['air_pollution_scaling'] / ((remind_air_pollution["ES|Transport|VKM|Pass|Road|LDV"] + remind_air_pollution["ES|Transport|VKM|Pass|Road|Bus"])*1000000000)


remind_air_pollution = remind_air_pollution['air_pollution_scaling']

# aggregate the air pollution to REMIND regions
air_pollution= air_pollution[air_pollution_variables]

#aggregation to the REMIND regions
region_mapping = pd.read_csv(region_mapping_path, sep=';')

air_pollution = air_pollution.merge(region_mapping[["CountryCode","RegionCode"]], left_on=['iso.iso.y'], right_on=['CountryCode'], how='left')

air_pollution = air_pollution.groupby(['RegionCode','year.1']).sum().reset_index()

#add a World row
air_pollution = pd.concat([air_pollution, air_pollution.groupby(['year.1']).mean().reset_index()])  
air_pollution['RegionCode'] = air_pollution['RegionCode'].replace(np.nan, "World")

#merge remind_air_pollution and air_pollution and scale with air pollution scaling

air_pollution = air_pollution.merge(remind_air_pollution, left_on=['RegionCode'], right_on=['Region'], how='left')
air_pollution = air_pollution.replace("World", 'GLO')

# scale with air pollution scaling
air_pollution['region'] = air_pollution["RegionCode"] 
air_pollution['period'] = air_pollution["year.1"] 

air_pollution['air_pollution_death'] = air_pollution["data_df_mort_pm.all_mort"] * air_pollution["air_pollution_scaling"]
air_pollution['air_pollution_DALY'] = air_pollution["daly_total"]* air_pollution["air_pollution_scaling"]
air_pollution['air_pollution_cost_indirect'] = air_pollution["indirect_daly_cost"]* air_pollution["air_pollution_scaling"]
air_pollution['air_pollution_cost_direct'] = air_pollution["direct_daly_cost"] *  air_pollution["air_pollution_scaling"]

air_pollution = air_pollution[['region','period','air_pollution_death','air_pollution_DALY','air_pollution_cost_indirect','air_pollution_cost_direct']]

# we scale the different technologies here with their LCA particualte matter shares and replace PMF
# we use a undifferentiated average across all technologies !!
# add a new row for air pollution deahts, dalys, indirect daly cost, direct daly cost

################ change this to period specific
data["average"] = data.mean(numeric_only=True, axis=1)

# calc the relative PMF compared to the average
data = data.append(data.filter(like = 'particulate', axis=0).rename(index={'particulate matter formation': 'air_pollution'})/  data.filter(like = 'particulate', axis=0)["average"].iloc[0]
)

# scale this factor to the air_pollution sclaing
data = data.transpose()
data = data.reset_index()


data["region"] = data['transport mode'].str.partition(', None')[0].str.strip().str[-0-3:]


#add year
data["period"] = data['transport mode'].str.partition(',')[0].str.partition("'")[2].str.partition("'")[0]

data['transport mode'] = data['transport mode'].str.partition(',')[2]

data['region'] = data['region'].replace('age', "GLO")
data['region'] = data['region'].replace(' CH', "EUR")

data['period'] = data['period'].replace('', "2050")
data['period'] = data['period'].astype(np.float64)

#left join data and air pollution

data = data.merge(air_pollution, on=['region','period'], how='left')

data['air_pollution_death'] = data['air_pollution_death'] * data['air_pollution']
data['air_pollution_DALY'] = data['air_pollution_DALY'] * data['air_pollution']
data['air_pollution_cost_indirect'] = data['air_pollution_cost_indirect'] * data['air_pollution']
data['air_pollution_cost_direct'] = data['air_pollution_cost_direct'] * data['air_pollution']
#drop scaling column
data = data.drop(columns=['air_pollution'])
#drop old LCA air pollution column
data = data.drop(columns=['particulate matter formation'])




