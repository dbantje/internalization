# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Calculate the technology specific LCA impact 
# # end point level
#
# 0. import brightway2 results for the coal exit scenarios
# 0. calcualte the energy market electricity high voltage mix RECIPE total, HH, ED and RD
# (indicator x scenarios x period)
#
# 0. calculate the mid to endpoint DALYS, and species.year
#
# 0. do the same for the other SE carriers 
#     1. coarse match of REMIND and brightway technologies on REMIND SE level
#     2. calculate the RECIPE total (3) for all these technologies in table form
#        (technologies x indicator x scenarios x period)
#
# 1. use these tables to do a post processing (R) to aggregate and monetize the DALYS and monetize the ED
#      1. read in results of this script
#      2. use conversion factors for ecoinvent to REMIND
#      3. use a REMIND activity mapping
#      4. monetize DALYS and PDFs
#     
# further developments
# 1. adjust mid to end pop with SSP pop data country specific?
# 2. adjust the land occupation indicator with country specific species data
# 3. disaggregate the land occupating indicator 
# 4. update to RECIPE 2016 when available
# 5. dont use the global mean for SE technologies
# 6. automatic scenario name mapping

# the se regional mapping is global

# # import the brightway2 results for the coal exit scenarios




import brightway2 as bw

from wurst.ecoinvent import filters
from constructive_geometries import Geomatcher

from wurst import searching as ws


from helpers import eimod, lcahelp, ei2rmnd

import pandas as pd
import numpy as np
from bw2data import Database





# take the 
#bw.restore_project_directory("data\\brightway\\brightway2-project-default-backup.09-April-2020-02-56PM.tar.gz")
bw.restore_project_directory("data\\brightway\\brightway2-project-default-backup.18-December-2020-10-06PM.tar.gz")

bw.projects.set_current("default")
project = 'default'

# ## create a scenarios dictionary

# +
database_dict = {}
database_dict['ecoinvent_added_technologies'] = {'year' : 2015 , 'scenario' : 'ecoinvent'}

for year in years:
    for scenario in scenarios:
        if year == 2015 and scenario != 'ecoinvent_added_technologies':
            continue
        db_name = 'ecoinvent_Remind_' + scenario + '_' + str(year)
        database_dict[db_name] = {'year' : year , 'scenario' : scenario}
# -

# ## calcualte the LCA for all scenarios and years for all SE technologies

# +
#select the modeled processes
#this is a crude mapping of of REMIND SE to ecoinvent progcesses
#we need to adjust the units later in the post processing to match EJ

# +
# see data\ecoinventtoREMINDtechnmap.csv for SE technology mapping
# we use the electricity generation processes as a base and check what the input is to adjust them accordingly
# see the A matrix of the elec processes for the input A.technology.final.technology.list.variable.csv
# ergo if a kWh hard coal needs 5kg of coal we:

# elec to raw input = /5
# 1kWh to 1EJ =  /3.6*10^12
# elec hard coal to se hard coal /3.6*10^12 /5kg

# this factor is in the conv column of the ecoinventtoREMINDtechmap.csv

# +
# se _filter filters for electricity high voltage and all the other SE ecoinvent processes mapped above except for transport

# +
############
# -

se_filter1=[ws.equals('name', 'heat and power co-generation, biogas, gas engine'),ws.equals('unit', 'megajoule')]
se_filter2=[ws.equals('name', 'heat and power co-generation, natural gas, conventional power plant, 100MW electrical'),ws.equals('unit', 'megajoule')]
se_filter3=[ws.equals('name', 'heat and power co-generation, hard coal'),ws.equals('unit', 'megajoule')]
se_filter4=[ws.equals('name', 'electricity production, deep geothermal'),ws.equals('unit', 'kilowatt hour')]
se_filter5=[ws.equals('name', 'hydrogen, liquid | market for hydrogen, liquid'),ws.equals('unit', 'kg')]
se_filter6=[ws.equals('name', 'heat and power co-generation, oil'),ws.equals('unit', 'megajoule')]
se_filter7=[ws.equals('name', 'heat production, softwood chips from forest, at furnace 300kW, state-of-the-art 2014'),ws.equals('unit', 'megajoule')]
se_filter8=[ws.equals('name', 'heat and power co-generation, biogas, gas engine'),ws.equals('unit', 'kilowatt hour')]
se_filter9=[ws.equals('name', 'electricity production, hard coal'),ws.equals('unit', 'kilowatt hour')]
se_filter10=[ws.equals('name', 'electricity production, natural gas, conventional power plant'),ws.equals('unit', 'kilowatt hour')]
se_filter11=[ws.equals('name', 'electricity production, oil'),ws.equals('unit', 'kilowatt hour')]
se_filter12=[ws.equals('name', 'electricity production, nuclear, pressure water reactor'),ws.equals('unit', 'kilowatt hour')]
se_filter13=[ws.equals('name', 'electricity production, hydro, run-of-river'),ws.equals('unit', 'kilowatt hour')]
se_filter14=[ws.equals('name', 'electricity production, solar thermal parabolic trough, 50 MW'),ws.equals('unit', 'kilowatt hour')]
se_filter15=[ws.equals('name', 'electricity production, photovoltaic, 570kWp open ground installation, multi-Si'),ws.equals('unit', 'kilowatt hour')]
se_filter16=[ws.equals('name', 'electricity production, wind, >3MW turbine, onshore'),ws.equals('unit', 'kilowatt hour')]

#### these processes only have global values in ecoinvent, we therfore treat them differently
ta_filter1=[ws.equals('name', 'market for transport, passenger car, medium size, diesel, EURO 5'),ws.equals('unit', 'kilometer')]
ta_filter2=[ws.equals('name', 'market for transport, passenger car, medium size, petrol, EURO 5'),ws.equals('unit', 'kilometer')]



electricity_markets = {}
se_technologies = {}
ta_technologies = {}


for db_name in database_dict:
    electricity_markets[db_name] = [x for x in ws.get_many(Database(db_name), *[*eimod.electricity_market_filter_high_voltage, ws.exclude(ws.contains('name', 'Swiss Federal Railways')),
                                                                                ws.exclude(ws.contains('name', 'label-certified'))])]																				
    se_technologies[db_name] = [x for x in ws.get_many(Database(db_name), *[*se_filter1])]+ [x for x in ws.get_many(Database(db_name), *[*se_filter2])] + [x for x in ws.get_many(Database(db_name), *[*se_filter3])]+[x for x in ws.get_many(Database(db_name), *[*se_filter4])]+[x for x in ws.get_many(Database(db_name), *[*se_filter5])]+[x for x in ws.get_many(Database(db_name), *[*se_filter6])]+[x for x in ws.get_many(Database(db_name), *[*se_filter7])]+[x for x in ws.get_many(Database(db_name), *[*se_filter8])]+[x for x in ws.get_many(Database(db_name), *[*se_filter9])]+[x for x in ws.get_many(Database(db_name), *[*se_filter10])]+[x for x in ws.get_many(Database(db_name), *[*se_filter11])]+[x for x in ws.get_many(Database(db_name), *[*se_filter12])]+[x for x in ws.get_many(Database(db_name), *[*se_filter13])]+[x for x in ws.get_many(Database(db_name), *[*se_filter14])]+[x for x in ws.get_many(Database(db_name), *[*se_filter15])]+[x for x in ws.get_many(Database(db_name), *[*se_filter16])]
    ta_technologies[db_name] = [x for x in ws.get_many(Database(db_name), *[*ta_filter1])] + [x for x in ws.get_many(Database(db_name), *[*ta_filter2])]


# -

# ## choose the RECIPE categories, we use 2008 with updated 2016 where direct matches are possible

# calculate all the ReCiPe Midpoint (H) 2008 version and multiply with mid to endpoints to get DALYS (HH), Species.year (ET) and maybe resources.

#method=[m for m in methods if "'ReCiPe Midpoint (H) V1.13'," in str(m)]
method = {
    'CC':('IPCC 2013', 'climate change', 'GWP 100a'),
    'TA':('ReCiPe Midpoint (H) V1.13', 'terrestrial acidification', 'TAP100'),
    'POF':('ReCiPe Midpoint (H) V1.13','photochemical oxidant formation','POFP'),
    'PMF':('ReCiPe Midpoint (H) V1.13', 'particulate matter formation', 'PMFP'),
    'MD':('ReCiPe Midpoint (H) V1.13', 'metal depletion', 'MDP'),
    'HT':('ReCiPe Midpoint (H) V1.13', 'human toxicity', 'HTPinf'),
    'MET':('ReCiPe Midpoint (H) V1.13', 'marine ecotoxicity', 'METPinf'),
    'ME':('ReCiPe Midpoint (H) V1.13', 'marine eutrophication', 'MEP'),
    'FD':('ReCiPe Midpoint (H) V1.13', 'fossil depletion', 'FDP'),
    'IR':('ReCiPe Midpoint (H) V1.13', 'ionising radiation', 'IRP_HE'),
    'OD':('ReCiPe Midpoint (H) V1.13', 'ozone depletion', 'ODPinf'),
    'FET':('ReCiPe Midpoint (H) V1.13', 'freshwater ecotoxicity', 'FETPinf'),
    'TET':('ReCiPe Midpoint (H) V1.13', 'terrestrial ecotoxicity', 'TETPinf'),
    'ALO':('ReCiPe Midpoint (H) V1.13', 'agricultural land occupation', 'ALOP'),
    'NLT':('ReCiPe Midpoint (H) V1.13', 'natural land transformation', 'NLTP'),
    'ULO':('ReCiPe Midpoint (H) V1.13', 'urban land occupation', 'ULOP'),
    'WD':('ReCiPe Midpoint (H) V1.13', 'water depletion', 'WDP'),
    'FE':('ReCiPe Midpoint (H) V1.13', 'freshwater eutrophication', 'FEP')}

# select the matching categories, we only use the IPCC CC instead of the RECIPE one, however they have the same unit so we can use the same mid to endpoint factor
#
# we use the RECIPE 1.11
#     1. we developed a land occupation method ourselves to include it in the Ecosystem toxicity
#     2. we did a crude matching of land occupation categories to urban and rural
#     3. we took the average of each of the two categories
#     4. for land transformation we use to and from forest and take the average terrestrial species per m2

categories = ['ALO', 'CC', 'FD', 'FET', 'FE',
              'HT', 'IR', 'MET', 'ME', 'MD',
              'NLT', 'OD', 'PMF', 'POF', 'TA',
              'TET', 'ULO', 'WD']
categories_total = ['ALO', 'CC', 'FD', 'FET', 'FE',
              'HT', 'IR', 'MET', 'ME', 'MD',
              'NLT', 'OD', 'PMF', 'POF', 'TA',
              'TET', 'ULO', 'WD',
                    'total'] 
					
categories = ['ALO']
categories_total = ['ALO','total'] 

# perform LCA for all elec markets

# +
electricity_markets_LCIA= {}
se_technologies_LCIA= {}
ta_technologies_LCIA= {}

for db_name in database_dict:
    for category in categories:
        electricity_markets_LCIA[(database_dict[db_name]['year'], database_dict[db_name]['scenario'], category)] = lcahelp.LCA_to_df(electricity_markets[db_name], cats = [category], names=['name', 'location', 'unit'])
		
"""		
        se_technologies_LCIA[(database_dict[db_name]['year'], database_dict[db_name]['scenario'], category)]  = lcahelp.LCA_to_df(se_technologies[db_name], cats = [category], names=['name', 'location', 'unit'])

for db_name in database_dict:
    for category in categories:
        ta_technologies_LCIA[(database_dict[db_name]['year'], database_dict[db_name]['scenario'], category)]  = lcahelp.LCA_to_df(ta_technologies[db_name], cats = [category], names=['name', 'location', 'unit'])

# 

# load the mid to end characterization factors
# calculate all results and aggregate them to DALYS and ET

# +
mid_end = pd.read_csv("data\\mapping\\recipe_mid_to_end.csv", sep = ';', index_col = 0) 


# -

end_categories = ['HH_H', 'ED_H', 'RA_H'] 
    

end_mid_electricity_markets = {}
end_mid_se_technologies = {}
end_mid_ta_technologies = {}
for db_name in database_dict:        
    for end_category in end_categories:
        for category in categories:
                    
            end_mid_electricity_markets[(database_dict[db_name]['year'], database_dict[db_name]['scenario'],  end_category, category)] = electricity_markets_LCIA[(database_dict[db_name]['year'], database_dict[db_name]['scenario'], category)] * mid_end.loc[category, end_category]
            end_mid_se_technologies[(database_dict[db_name]['year'], database_dict[db_name]['scenario'],  end_category, category)]     = se_technologies_LCIA[(database_dict[db_name]['year'], database_dict[db_name]['scenario'], category)] * mid_end.loc[category, end_category]  
            end_mid_ta_technologies[(database_dict[db_name]['year'], database_dict[db_name]['scenario'],  end_category, category)]     = ta_technologies_LCIA[(database_dict[db_name]['year'], database_dict[db_name]['scenario'], category)] * mid_end.loc[category, end_category]  			



end_end_electricity_markets = pd.concat(end_mid_electricity_markets, names=['period','scenario','end_category', 'category', 'variable'])
end_end_se_technologies     = pd.concat(end_mid_se_technologies, names=['period','scenario','end_category', 'category', 'variable'])
end_end_ta_technologies     = pd.concat(end_mid_ta_technologies, names=['period','scenario','end_category', 'category', 'variable'])

# sum up all end categories

end_end_electricity_markets = end_end_electricity_markets.groupby(['period','scenario', 'end_category', 'variable']).sum()
end_end_se_technologies     = end_end_se_technologies.groupby(['period','scenario', 'end_category', 'variable']).sum()
end_end_ta_technologies     = end_end_ta_technologies.groupby(['period','scenario', 'end_category', 'variable']).sum()

end_end_electricity_markets['total'] = end_end_electricity_markets.sum(axis = 1)
end_end_se_technologies['total']     = end_end_se_technologies.sum(axis = 1)
end_end_ta_technologies['total']     = end_end_ta_technologies.sum(axis = 1)

end_end_electricity_markets.to_csv('data\\interm\\lca_electricity_markets_per_tec' + project + '.csv', sep=';')
end_end_se_technologies.to_csv('data\\interm\\lca_se_technologies_per_tec' + project + '.csv', sep=';')
end_end_ta_technologies.to_csv('data\\interm\\lca_ta_technologies_per_tec' + project + '.csv', sep=';')

end_end_electricity_markets  = pd.read_csv('data\\interm\\lca_electricity_markets_per_tec' + project + '.csv', sep = ';') 
end_end_se_technologies      = pd.read_csv('data\\interm\\lca_se_technologies_per_tec' + project + '.csv', sep = ';') 
end_end_ta_technologies      = pd.read_csv('data\\interm\\lca_ta_technologies_per_tec' + project + '.csv', sep = ';') 

# add location

end_end_electricity_markets['ecoinvent'] = end_end_electricity_markets['variable'].str.split("'").str[3].str.split("-").str[0]
end_end_se_technologies['ecoinvent']     = end_end_se_technologies['variable'].str.split("'").str[3].str.split("-").str[0]
end_end_ta_technologies['ecoinvent']     = end_end_ta_technologies['variable'].str.split("'").str[3].str.split("-").str[0]


end_end_electricity_markets['variable'].str.split(",").str[0].str.split("'").str[1] +", " + end_end_electricity_markets['variable'].str.split("'").str[5]

# match REMIND locations

end_end_electricity_markets['variable'] = end_end_electricity_markets['variable'].str.split(",").str[0].str.split("'").str[1] +", " + end_end_electricity_markets['variable'].str.split("'").str[5]
end_end_se_technologies['variable'] = end_end_se_technologies['variable'].str.split("'").str[1]+", " + end_end_se_technologies['variable'].str.split("'").str[5]
end_end_ta_technologies['variable'] = end_end_ta_technologies['variable'].str.split("'").str[1]+", " + end_end_ta_technologies['variable'].str.split("'").str[5]


region_mapping_elec     = pd.read_csv("data\\mapping\\regionmappingREMIND_noFullISO_elec.csv", sep = ';', index_col = 0)
region_mapping_se_world = pd.read_csv("data\\mapping\\regionmappingREMIND_noFullISO_elec_world.csv", sep = ';', index_col = 0)
region_mapping_ta_world = pd.read_csv("data\\mapping\\regionmappingREMIND_noFullISO_world.csv", sep = ';', index_col = 0)


end_end_electricity_markets = pd.merge(end_end_electricity_markets, region_mapping_elec, on= "ecoinvent", how= "inner")
end_end_se_technologies_world = pd.merge(end_end_se_technologies, region_mapping_se_world, on= "ecoinvent", how= "inner")
end_end_se_technologies = pd.merge(end_end_se_technologies, region_mapping_elec, on= "ecoinvent", how= "inner")


end_end_ta_technologies_world = pd.merge(end_end_ta_technologies, region_mapping_se_world, on= "ecoinvent", how= "inner")
end_end_ta_technologies = pd.merge(end_end_ta_technologies, region_mapping_ta_world, on= "ecoinvent", how= "inner")

# average over multiple ecoinvent locations, we have India with several grids for example


end_end_electricity_markets = end_end_electricity_markets.groupby(['period','scenario', 'end_category', 'variable', 'RegionCode']).mean()
end_end_se_technologies = end_end_se_technologies.groupby(['period','scenario', 'end_category', 'variable', 'RegionCode']).mean()
end_end_se_technologies_world = end_end_se_technologies_world.groupby(['period','scenario', 'end_category', 'variable', 'RegionCode']).mean()

end_end_ta_technologies = end_end_ta_technologies.groupby(['period','scenario', 'end_category', 'variable', 'RegionCode']).mean()
end_end_ta_technologies_world = end_end_ta_technologies_world.groupby(['period','scenario', 'end_category', 'variable', 'RegionCode']).mean()

#
#

end_end_se_technologies.to_csv('data\\interm\\lca.csv', sep=';')
end_end_se_technologies = pd.read_csv("data\\interm\\lca.csv", sep = ';') 

end_end_se_technologies_world.to_csv('data\\interm\\lca.csv', sep=';')
end_end_se_technologies_world = pd.read_csv("data\\interm\\lca.csv", sep = ';') 

end_end_ta_technologies.to_csv('data\\interm\\lca.csv', sep=';')
end_end_ta_technologies = pd.read_csv("data\\interm\\lca.csv", sep = ';') 
end_end_ta_technologies_world.to_csv('data\\interm\\lca.csv', sep=';')
end_end_ta_technologies_world = pd.read_csv("data\\interm\\lca.csv", sep = ';') 

end_end_electricity_markets.to_csv('data\\interm\\lca.csv', sep=';')
end_end_electricity_markets = pd.read_csv("data\\interm\\lca.csv", sep = ';')

# +
#end_end_se_technologies_world = pd.merge(end_end_se_technologies, region_mapping_elec,  left_on=['RegionCode'], right_on=['World'], how= "right")

# +
#end_end_se_technologies['RegionCode'] = end_end_se_technologies['Remind']
#delte World and -- column
#end_end_se_technologies = end_end_se_technologies.drop('World', 1)
#end_end_se_technologies = end_end_se_technologies.drop('World', 1)
#end_end_se_technologies = end_end_se_technologies.drop('Remind', 1)
# -

end_end_se_technologies = end_end_se_technologies.append(end_end_se_technologies_world)
end_end_ta_technologies = end_end_ta_technologies.append(end_end_ta_technologies_world)

end_end = end_end_electricity_markets.append(end_end_se_technologies.append(end_end_ta_technologies))

pd.DataFrame(end_end)

# join with REMIND mapping and energy content conversion

# +
# join with REMIND mapping
## read REMIND mapping
## join
## divide with convertion to go to EJ
## the conversion is the result of the ecoinvent process unit conversion, the SE per process input and and assumed engine efficiency for the process itself.

# on the example of biogas:
#ecoinvent	ecoinvent_unit	REMIND	REMIND_unit	conv_unit	SE_eco_input	engine_eff	conv
#heat and power co-generation, biogas, gas engine	kWh	SE|Biomass	EJ	3.60E-12	0.344325	0.5	5.23E-12

#kWh to EJ, SE biogas input in process, and engine_eff for SE to heat or similar out of biogas.
# is the engine efficiency necessary here or is it double counted from the REMIND efficiency factors?

# -

tech_mapping = pd.read_csv("data\\mapping\\ecoinventtoREMINDtechmap_new.csv", sep = ';', index_col = False) 

end_end = pd.merge(end_end, tech_mapping, left_on=['variable'], right_on=['ecoinvent'], how= "left")

end_end.to_csv('results\\lca_per_tec_prelim_' + project  + '.csv', sep=';')

end_end['total'] = end_end['total'] / end_end['conv']
end_end[categories]  = end_end[categories].divide(end_end["conv"], axis="index")

# cleanup

end_end_out = end_end[['period','scenario','end_category', 'RegionCode', 'variable']]
end_end_out[categories] = end_end[categories]
end_end_out['total'] = end_end['total']

# +
#end_end_out = pd.melt(end_end_out, id_vars=['period', 'scenario','RegionCode','end_category','variable'], value_vars=categories_total)
# -

# write results

#end_end_out.to_csv('results\\lca_per_tec_' + project  + '.csv', sep=';')
#end_end_out = pd.read_csv("results\\lca_per_tec_default.csv", sep = ';') 

# # read all REMIND results and calculate the actual model specific LCA including the monetization
# 1. read the vm_prodse produced from the vm_prodse_gdx.R script
# 2. join with end_end_out and multiply
# 3. join with vsl_coeff
#     1. join vsl_coeff with regions mapping and average over region

# rbind se and fe

vm_prod_se = pd.read_csv("data\\Remind output files\\vm_prodse_internalization.csv", sep = ',', index_col = 0) 
vm_prod_fe = pd.read_csv("data\\Remind output files\\vm_prodfe_internalization.csv", sep = ',', index_col = 0) 

vm_prod_fe.rename(columns={'ttot': 'tall'}, inplace=True)

vm_prod = pd.concat([vm_prod_se, vm_prod_fe])

# define colnames that are years

# rename columns

end_end_out=end_end_out.melt(id_vars=['period','scenario','end_category','RegionCode', 'variable'],var_name='category')

vm_prod.columns=['scenario','period','region','tech','SE']
end_end_out.columns=['period','scenario','end_category','region','tech','category','lca']

#add a world row with the total
vm_prod_world = pd.pivot_table(vm_prod, values='SE', index=['scenario','period','tech'], aggfunc=np.sum, fill_value=0)
vm_prod_world.reset_index(level=['scenario', 'period', 'tech'], inplace=True)
vm_prod_world
vm_prod_world['region'] = 'World'
vm_prod = vm_prod.append(vm_prod_world)

# join LCA results with the vm_prodse

# map the technologies from ecoinvent to vm_prodSe

vm_prodse_ecoinvent_tech_mapping     = pd.read_csv("data\\mapping\\REMIND_reportingtechmap_new.csv", sep = ';', index_col = 0)

vm_prodse_ecoinvent_tech_mapping = vm_prodse_ecoinvent_tech_mapping.dropna()

vm_prod = pd.merge(vm_prod, vm_prodse_ecoinvent_tech_mapping , left_on=['tech'],right_on=['REMIND'], how= "left")

vm_prod['period']=vm_prod['period'].astype(int)
end_end_out['period']=end_end_out['period'].astype(int)



# map scenarios

end_end_out['scenario'][end_end_out['scenario'] == 'ecoinvent'] = 'Base'

end_end_out.rename(columns = {'tech':'ecoinvent'}, inplace = True) 

#I guess the error is here, check if its scenario mapping or what?


end_end_out = pd.merge(vm_prod, end_end_out , left_on=['period','scenario','ecoinvent','region'],right_on=['period','scenario','ecoinvent','region'], how= "left")

#filter for years
end_end_out = end_end_out[end_end_out['period'].isin(years)]



# +
# if lca is missing for region, tech, period take the world values
# need to write all impact categories
# filter out the zero columns
# filter the world and merge the two
# append this to the original one


end_end_out_missing = end_end_out[end_end_out['lca'].isna()]
end_end_out_world  = end_end_out[end_end_out['region'] == 'World']
#end_end_out_missing.to_csv('test.csv')
#end_end_out_missing

# +
# delete same columns from missing
# delete SE from world and "world" from world

end_end_out_missing = pd.merge(end_end_out_missing,end_end_out_world , left_on=['scenario','period','tech'],right_on=['scenario','period','tech'], how= "right")##


# -

end_end_out_missing.drop(columns=['region_y',
'SE_y',
'descripton_x',
'ecoinvent_process_x',
'ecoinvent_unit_x',
'ecoinvent_x',
'conv_unit_x',
'SE_eco_input_x',
'conv_x',
'end_category_x',
'category_x',
'lca_x'], axis = 1, inplace = True)

# +
end_end_out_missing.rename(columns={     'region_x':'region','SE_x':'SE', 'descripton_y':'descripton',
                                         'ecoinvent_process_y':'ecoinvent_process', 'ecoinvent_unit_y':'ecoinvent_unit',
                                          'ecoinvent_y':'ecoinvent', 'conv_unit_y':'conv_unit',
                                          'SE_eco_input_y':'SE_eco_input','conv_y':'conv_y',
                                          'end_category_y':'end_category','category_y':'category',
                                          'lca_y':'lca'}, inplace = True)#



# +
end_end_out = end_end_out[end_end_out['lca'].notna()].append(end_end_out_missing)


# -

end_end_out.drop(columns=['descripton',
'ecoinvent_process',
'ecoinvent_unit',
'ecoinvent',
'conv_y',
'conv_unit'],axis = 1, inplace = True)

end_end_out['SE'] = end_end_out['SE'].replace(0,10**-9)
end_end_out['lca_impact'] = end_end_out['lca'] * end_end_out['SE']

project = 'project'

# write results

end_end_out.to_csv('results\\lca_REMIND_' + project  + '.csv', sep=';')

# ## read the VSL and multiply the HH and ED values with the VSL and end category

# read the vsl coeff and do the average over REMIND regions
#     map to REMIND regions
# join vsl with end_end_out 
# multiply with monetization for ED and HH and R?

vsl_coeff = pd.read_csv("data//vsl_coeff.csv", sep = ';', index_col = 0) 

# join with REMIND region mapping

# ## condense to REMIND region results

vsl_coeff = pd.merge(vsl_coeff, region_mapping_elec, left_on= "Region", right_on='CountryCode', how= "left")
vsl_coeff = vsl_coeff.groupby(['RegionCode','Year','Data2']).mean()
vsl_coeff.columns = ['vsl_coeff']

# ## join end_end_out and vsl_coeff

end_end_out = pd.merge(end_end_out, vsl_coeff, left_on= ["period","region"], right_on=["Year","RegionCode"], how= "left")

# ## add the monetization

# #VSL from OECD (2012) Mortality risk valuation in environment, health and transport policies. OECD Publishing
# #divded by the 2017 ration of mortality to DALYS, all cause all ages
# u_DALY <- 3600000 / 25 [$2005]

# +
#HH_H [DALYS] = 144000  

#ED_H [species.yr] = 50447778 , Kiuk 2008 gives a value from the NEEDs project of 0.25 €2004/PDF/m2/year for the EU 25 based on restauration cost of from Germany to go from integreated arable to organic arable (the lowest cost to transform land)

#@Ott2007 gives a spatially weighted average to transform "Built Up Land" into
#other land for EU-25: 1.625 $2005/PDF/m2

#this equals 0.0165 $2005/PDF/m2/year if we apply 100a for land-use change to last

#terrestrial species density: 1.48 E-8 [1/m2][3] global data from recipe
#	CF PDF = SD PDF, only terrestrial
#	CF [↕species.yr] = PDF.m2.yr * species/m²
#	€/ CF [species.yr]= CF [↕species.yr] *1/( species/m²)
#1/1.48 E-8 [1/m2] = 67567568 [1/(species/m²)]
#0.0165$ [2005] per m² per year * 67567568 [1/(species/m²)] = 1114864.872 $/species year


#RD_H [$] = 1.09

# +
end_end_out['HH_H_marginal_cost'] = 1
end_end_out['ED_H_marginal_cost'] = 11148648.72                                    
end_end_out['RA_H_marginal_cost'] = 1.09

end_end_out['HH_H_cost_total']=0
end_end_out['ED_H_cost_total']=0
end_end_out['RA_H_cost_total']=0
# -



# +
end_end_out['HH_H_cost_total'][end_end_out['end_category'] == 'HH_H'] = end_end_out['HH_H_marginal_cost'] * end_end_out['lca_impact']
end_end_out['ED_H_cost_total'][end_end_out['end_category'] == 'ED_H'] = end_end_out['ED_H_marginal_cost'] * end_end_out['lca_impact']
end_end_out['RA_H_cost_total'][end_end_out['end_category'] == 'RA_H'] = end_end_out['RA_H_marginal_cost'] * end_end_out['lca_impact']

end_end_out['HH_H_cost_total_vsl']=0
end_end_out['ED_H_cost_total_vsl']=0
end_end_out['RA_H_cost_total_vsl']=0
end_end_out['HH_H_cost_total_vsl'][end_end_out['end_category'] == 'HH_H'] = end_end_out['HH_H_cost_total'] * end_end_out['vsl_coeff']
end_end_out['ED_H_cost_total_vsl'][end_end_out['end_category'] == 'ED_H'] = end_end_out['ED_H_cost_total'] * end_end_out['vsl_coeff']
end_end_out['RA_H_cost_total_vsl'][end_end_out['end_category'] == 'RA_H'] = end_end_out['RA_H_cost_total'] * end_end_out['vsl_coeff']
# -

# ## write all categories and total

end_end_out.to_csv('results\\lca_REMIND_monetization_' + project  + '.csv', sep=';')

# +
#run the policy_cost.R to get the mitigation cost

# run the calclca_cost_tax.R to add the AP to the LCA and calucalte tax, here the mitigation cost are added
# -
"""