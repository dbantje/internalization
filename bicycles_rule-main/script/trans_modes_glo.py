
#read remind results
remind = pd.read_csv(remind_report, sep=';')
remind = remind[['Model', 'Scenario', 'Region', 'Variable', 'Unit',str(year)]]

#extract all transport modes and the global
remind_pass_road_ldv_total                   = remind.query('Variable == ("ES|Transport|Pass|Road|LDV") and Region.str.contains("World")', engine='python')
remind_pass_road_ldv_bev_total               = remind.query('Variable == ("ES|Transport|Pass|Road|LDV|BEV") and Region.str.contains("World")', engine='python')
remind_pass_road_ldv_fcev_total              = remind.query('Variable == ("ES|Transport|Pass|Road|LDV|FCEV") and Region.str.contains("World")', engine='python')
remind_pass_road_ldv_gases_total             = remind.query('Variable == ("ES|Transport|Pass|Road|LDV|Gases") and Region.str.contains("World")', engine='python')
remind_pass_road_ldv_hybrid_electric_total   = remind.query('Variable == ("ES|Transport|Pass|Road|LDV|Hybrid Electric") and Region.str.contains("World")', engine='python')
remind_pass_road_ldv_hybrid_liquids_total    = remind.query('Variable == ("ES|Transport|Pass|Road|LDV|Hybrid Liquids") and Region.str.contains("World")', engine='python')
remind_pass_road_ldv_liquids_total           = remind.query('Variable == ("ES|Transport|Pass|Road|LDV|Liquids") and Region.str.contains("World")', engine='python')

#select single processes
remind_pass_road_ldv_total             = remind.query('Variable == ("ES|Transport|Pass|Road|LDV") and Region.str.contains("World")', engine='python')
remind_pass_road_ldv_bev               = remind.query('Variable == ("ES|Transport|Pass|Road|LDV|BEV") and Region != ("World")', engine='python')
remind_pass_road_ldv_fcev              = remind.query('Variable == ("ES|Transport|Pass|Road|LDV|FCEV") and Region != ("World")', engine='python')
remind_pass_road_ldv_gases             = remind.query('Variable == ("ES|Transport|Pass|Road|LDV|Gases") and Region != ("World")', engine='python')
remind_pass_road_ldv_hybrid_electric   = remind.query('Variable == ("ES|Transport|Pass|Road|LDV|Hybrid Electric") and Region != ("World")', engine='python')
remind_pass_road_ldv_hybrid_liquids    = remind.query('Variable == ("ES|Transport|Pass|Road|LDV|Hybrid Liquids") and Region != ("World")', engine='python')
remind_pass_road_ldv_liquids           = remind.query('Variable == ("ES|Transport|Pass|Road|LDV|Liquids") and Region != ("World")', engine='python')

#construct a GLO total:
remind_pass_road_ldv_total = pd.concat([remind_pass_road_ldv_bev_total,remind_pass_road_ldv_fcev_total,remind_pass_road_ldv_gases_total,remind_pass_road_ldv_hybrid_electric_total, remind_pass_road_ldv_hybrid_liquids_total,remind_pass_road_ldv_liquids_total])

remind_pass_road_ldv_total['Variable'] = remind_pass_road_ldv_total['Variable'].map(ecoinvent_activity_transport_mapping)




#calculate regional share of global for each transport mode
remind_pass_road_ldv_bev[str(year)] = remind_pass_road_ldv_bev[str(year)].div(remind_pass_road_ldv_bev_total[str(year)].iloc[0])
remind_pass_road_ldv_fcev[str(year)] = remind_pass_road_ldv_fcev[str(year)].div(remind_pass_road_ldv_fcev_total[str(year)].iloc[0])
remind_pass_road_ldv_gases[str(year)] = remind_pass_road_ldv_gases[str(year)].div(remind_pass_road_ldv_gases_total[str(year)].iloc[0])
remind_pass_road_ldv_hybrid_electric[str(year)] = remind_pass_road_ldv_hybrid_electric[str(year)].div(remind_pass_road_ldv_hybrid_electric_total[str(year)].iloc[0])
remind_pass_road_ldv_hybrid_liquids[str(year)] = remind_pass_road_ldv_hybrid_liquids[str(year)].div(remind_pass_road_ldv_hybrid_liquids_total[str(year)].iloc[0])
remind_pass_road_ldv_liquids[str(year)] = remind_pass_road_ldv_liquids[str(year)].div(remind_pass_road_ldv_liquids_total[str(year)].iloc[0])


remind_pass_road_ldv_bev['Variable'] = remind_pass_road_ldv_bev['Variable'].map(ecoinvent_activity_transport_mapping)
remind_pass_road_ldv_fcev['Variable'] = remind_pass_road_ldv_fcev['Variable'].map(ecoinvent_activity_transport_mapping)
remind_pass_road_ldv_gases['Variable'] = remind_pass_road_ldv_gases['Variable'].map(ecoinvent_activity_transport_mapping)
remind_pass_road_ldv_hybrid_electric['Variable'] = remind_pass_road_ldv_hybrid_electric['Variable'].map(ecoinvent_activity_transport_mapping)
remind_pass_road_ldv_hybrid_liquids['Variable'] = remind_pass_road_ldv_hybrid_liquids['Variable'].map(ecoinvent_activity_transport_mapping)
remind_pass_road_ldv_liquids['Variable'] = remind_pass_road_ldv_liquids['Variable'].map(ecoinvent_activity_transport_mapping)

#replace region with the remind['Region']

remind_pass_road_ldv_bev['Variable']= [x.replace('region', str(y)) for x, y  in remind_pass_road_ldv_bev[['Variable','Region']].to_numpy()]
remind_pass_road_ldv_fcev['Variable']= [x.replace('region', str(y)) for x, y  in remind_pass_road_ldv_fcev[['Variable','Region']].to_numpy()]
remind_pass_road_ldv_gases['Variable']= [x.replace('region', str(y)) for x, y  in remind_pass_road_ldv_gases[['Variable','Region']].to_numpy()]
remind_pass_road_ldv_hybrid_electric['Variable']= [x.replace('region', str(y)) for x, y  in remind_pass_road_ldv_hybrid_electric[['Variable','Region']].to_numpy()]
remind_pass_road_ldv_hybrid_liquids['Variable']= [x.replace('region', str(y)) for x, y  in remind_pass_road_ldv_hybrid_liquids[['Variable','Region']].to_numpy()]
remind_pass_road_ldv_liquids['Variable']= [x.replace('region', str(y)) for x, y  in remind_pass_road_ldv_liquids[['Variable','Region']].to_numpy()]







# construct GLO activities from the regional share############################
# not very elegant, make a loop...
activities=[]

for index, row in remind_pass_road_ldv_bev.iterrows():
    variable = row['Variable'].split("'",2)[1]
    region=row['Region']
    activities.append([act for act in lca_databases[str(year)] if act['name']==str(variable) and act['location'] == str(region)][0])

#create a new World activity for every transport mode 
pass_road_ldv_bev = lca_databases[str(year)].new_activity(code = "pass_road_ldv_bev", name = str("pass_road_ldv_bev"), unit = "vkm")

#loop over regions
for i in range(len(activities)-1):
    activity = str(activities[i])
    print(activity)

    amt = remind_pass_road_ldv_bev.query("Variable == @activity")[str(year)]
    amt = pd.to_numeric(amt, errors='coerce')
    region= remind_pass_road_ldv_bev.query("Variable == @activity")['Region']
    print(amt.values[0])

    pass_road_ldv_bev.new_exchange(input=activities[i].key,amount=amt.values[0],unit='vkm', location=str(region),type='technosphere').save()

pass_road_ldv_bev.save()
       

activities=[]

for index, row in remind_pass_road_ldv_fcev.iterrows():
    variable = row['Variable'].split("'",2)[1]
    region=row['Region']
    activities.append([act for act in lca_databases[str(year)] if act['name']==str(variable) and act['location'] == str(region)][0])

#create a new World activity for every transport mode 
pass_road_ldv_fcev = lca_databases[str(year)].new_activity(code = "pass_road_ldv_fcev", name = str("pass_road_ldv_fcev"), unit = "vkm")

#loop over regions
for i in range(len(activities)-1):
    activity = str(activities[i])
    print(activity)

    amt = remind_pass_road_ldv_fcev.query("Variable == @activity")[str(year)]
    amt = pd.to_numeric(amt, errors='coerce')
    region= remind_pass_road_ldv_fcev.query("Variable == @activity")['Region']
    print(amt.values[0])

    pass_road_ldv_fcev.new_exchange(input=activities[i].key,amount=amt.values[0],unit='vkm', location=str(region),type='technosphere').save()

pass_road_ldv_fcev.save()


activities=[]

for index, row in remind_pass_road_ldv_gases.iterrows():
    variable = row['Variable'].split("'",2)[1]
    region=row['Region']
    activities.append([act for act in lca_databases[str(year)] if act['name']==str(variable) and act['location'] == str(region)][0])

#create a new World activity for every transport mode 
pass_road_ldv_gases = lca_databases[str(year)].new_activity(code = "pass_road_ldv_gases", name = str("pass_road_ldv_gases"), unit = "vkm")

#loop over regions
for i in range(len(activities)-1):
    activity = str(activities[i])
    print(activity)

    amt = remind_pass_road_ldv_gases.query("Variable == @activity")[str(year)]
    amt = pd.to_numeric(amt, errors='coerce')
    region= remind_pass_road_ldv_gases.query("Variable == @activity")['Region']
    print(amt.values[0])

    pass_road_ldv_gases.new_exchange(input=activities[i].key,amount=amt.values[0],unit='vkm', location=str(region),type='technosphere').save()

pass_road_ldv_gases.save()

activities=[]

for index, row in remind_pass_road_ldv_hybrid_electric.iterrows():
    variable = row['Variable'].split("'",2)[1]
    region=row['Region']
    activities.append([act for act in lca_databases[str(year)] if act['name']==str(variable) and act['location'] == str(region)][0])

#create a new World activity for every transport mode 
pass_road_ldv_hybrid_electric = lca_databases[str(year)].new_activity(code = "pass_road_ldv_hybrid_electric", name = str("pass_road_ldv_hybrid_electric"), unit = "vkm")

#loop over regions
for i in range(len(activities)-1):
    activity = str(activities[i])
    print(activity)

    amt = remind_pass_road_ldv_hybrid_electric.query("Variable == @activity")[str(year)]
    amt = pd.to_numeric(amt, errors='coerce')
    region= remind_pass_road_ldv_hybrid_electric.query("Variable == @activity")['Region']
    print(amt.values[0])

    pass_road_ldv_hybrid_electric.new_exchange(input=activities[i].key,amount=amt.values[0],unit='vkm', location=str(region),type='technosphere').save()

pass_road_ldv_hybrid_electric.save()

activities=[]

for index, row in remind_pass_road_ldv_hybrid_liquids.iterrows():
    variable = row['Variable'].split("'",2)[1]
    region=row['Region']
    activities.append([act for act in lca_databases[str(year)] if act['name']==str(variable) and act['location'] == str(region)][0])

#create a new World activity for every transport mode 
pass_road_ldv_hybrid_liquids = lca_databases[str(year)].new_activity(code = "pass_road_ldv_hybrid_liquids", name = str("pass_road_ldv_hybrid_liquids"), unit = "vkm")

#loop over regions
for i in range(len(activities)-1):
    activity = str(activities[i])
    print(activity)

    amt = remind_pass_road_ldv_hybrid_liquids.query("Variable == @activity")[str(year)]
    amt = pd.to_numeric(amt, errors='coerce')
    region= remind_pass_road_ldv_hybrid_liquids.query("Variable == @activity")['Region']
    print(amt.values[0])

    pass_road_ldv_hybrid_liquids.new_exchange(input=activities[i].key,amount=amt.values[0],unit='vkm', location=str(region),type='technosphere').save()

pass_road_ldv_hybrid_liquids.save()

activities=[]

for index, row in remind_pass_road_ldv_liquids.iterrows():
    variable = row['Variable'].split("'",2)[1]
    region=row['Region']
    activities.append([act for act in lca_databases[str(year)] if act['name']==str(variable) and act['location'] == str(region)][0])

#create a new World activity for every transport mode 
pass_road_ldv_liquids = lca_databases[str(year)].new_activity(code = "pass_road_ldv_liquids", name = str("pass_road_ldv_liquids"), unit = "vkm")

#loop over regions
for i in range(len(activities)-1):
    activity = str(activities[i])
    print(activity)

    amt = remind_pass_road_ldv_liquids.query("Variable == @activity")[str(year)]
    amt = pd.to_numeric(amt, errors='coerce')
    region= remind_pass_road_ldv_liquids.query("Variable == @activity")['Region']
    print(amt.values[0])

    pass_road_ldv_liquids.new_exchange(input=activities[i].key,amount=amt.values[0],unit='vkm', location=str(region),type='technosphere').save()

pass_road_ldv_liquids.save()
# construct GLO activities from the regional share############################

