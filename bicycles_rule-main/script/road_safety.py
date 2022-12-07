#Road safety, find justification why to assign deaths to cars
# seems obvious if you think about it, there are probably very few deaths of cyclist - pedestrian or cyclist - cyclist
# is it sensible to assign all road injuries to cars, maybe only pick car injuries?

techs = ["Cyclist road injuries",
        "Motor vehicle road injuries",
        "Motorcyclist road injuries",
        "Pedestrian road injuries"]


df_acc            =  df.loc[(df['rei_name'] == 'All risk factors')& (df['measure_name'] == 'DALY') & (df['cause_name'].isin(techs))]


# we need a key of how many bicycle and pedestrian injuries are caused by cars
# numbers are hard to find, Berlin police reported that in 2020 50% of bicycle accidents where not caused by the cyclists
#https://www.berlin.de/polizei/aufgaben/verkehrssicherheit/verkehrsunfallstatistik/

# nsc: US, of the 1,089 bicyclist deaths in 2019, 712 died in motor-vehicle crashes and 377 in other crashes
# this does not talk about whos fault it is but not sure if this is relevant since the car caused the deaths,
# germans say in 75% of all crashes that involve cars, it is mainly the cars fault
#https://www.sazbike.de/markt-politik/adfc-allgemeiner-deutscher/adfc-fahrradunfaelle-nehmen-drastisch-zu-1728531.html

# no car no deaths, maybe only a fall of the cyclist -> all crashes involing cars are assigned to cars
# pedestrians are assigned to cars

# 

edge_internalization_mapping_acc             = {"Cycle": "Cyclist road injuries",
                                                "Compact Car": "Motor vehicle road injuries",
                                               'Large Car': "Motor vehicle road injuries",
                                               'Large Car and SUV': "Motor vehicle road injuries",
                                               'Light Truck and SUV': "Motor vehicle road injuries",
                                               'Midsize Car': "Motor vehicle road injuries",
                                               'Mini Car': "Motor vehicle road injuries",
                                               'Moped': 'Motorcyclist road injuries',
                                               'Motorcycle (50-250cc)': 'Motorcyclist road injuries',
                                               'Motorcycle (>250cc)': 'Motorcyclist road injuries',
                                               'Subcompact Car': "Motor vehicle road injuries",
                                               'Van': "Motor vehicle road injuries",
                                               'Bus_tmp_vehicletype': "Motor vehicle road injuries",
                                               'Truck (0-3.5t)': "Motor vehicle road injuries",
                                               'Truck (18t)': "Motor vehicle road injuries",
                                               'Truck (26t)': "Motor vehicle road injuries",
                                               'Truck (40t)': "Motor vehicle road injuries",
                                               'Truck (7.5t)': "Motor vehicle road injuries"}

  

df_acc = df_acc.pivot(index=['year','remind','measure_name','rei_name'],columns = 'cause_name', values="value")

df_acc['Motor vehicle road injuries'] = df_acc['Motor vehicle road injuries'] + df_acc['Cyclist road injuries']* 0.654
df_acc['Motor vehicle road injuries'] = df_acc['Motor vehicle road injuries'] + df_acc['Pedestrian road injuries']
df_acc['Cyclist road injuries']       = df_acc['Cyclist road injuries']* 0.346
df_acc['Pedestrian road injuries']    = df_acc['Pedestrian road injuries']* 0

df_acc=df_acc.reset_index()

df_acc = pd.melt(df_acc, id_vars=['year','remind','rei_name','measure_name'],
        var_name="cause_name", 
        value_name="value",
        value_vars=['Cyclist road injuries', 'Motor vehicle road injuries','Motorcyclist road injuries','Pedestrian road injuries'])
		
REMIND = ["ES|Transport|VKM|Pass|Road|LDV",
          "ES|Transport|Pass|Road|Non-Motorized",
          "ES|Transport|VKM|Pass|Road|LDV|Two-Wheelers"]


#pedestrian is assigned to cars, its therefore missing in the mapping
GDB_REMIND_mapping                          = {"ES|Transport|Pass|Road|Non-Motorized": "Cyclist road injuries",
                                               "ES|Transport|VKM|Pass|Road|LDV": "Motor vehicle road injuries",
                                               "ES|Transport|VKM|Pass|Road|LDV|Two-Wheelers": "Motorcyclist road injuries"}
  
#filter the remind results for the activities abvove
remind_gbd       = remind.loc[(remind['Variable'].isin(REMIND))]

#map
remind_gbd['Variable'] = remind_gbd['Variable'].map(GDB_REMIND_mapping)
remind_gbd = pd.melt(remind_gbd, id_vars=['Model','Scenario','Region' ,'Variable','Unit'], var_name='Period', value_name='value',ignore_index=True)

#delete 
remind_gbd    = remind_gbd.loc[(~remind_gbd['Period'].isin(['Unnamed: 24']))]
remind_gbd['Period'] = remind_gbd['Period'].astype(int)


#merge
df_acc_remind = df_acc.merge(remind_gbd,  how='left',left_on = ['remind','year','cause_name'],right_on = ['Region','Period','Variable'])

df_acc_remind['value'] = df_acc_remind['value_x'] / df_acc_remind['value_y']
df_acc_remind['Unit']  = 'Road injuries per capita DALY/bn vkm'
df_acc_remind = df_acc_remind.dropna()
df_acc_remind = df_acc_remind[['Region','year','Scenario','Variable','Unit','value']]
df_acc_remind = df_acc_remind.rename(columns={"year":"Period"})
df_acc_remind
		