#low physical activity

# assigne confidence levels?

#Physical activity:

## mostly rather old affected, maybe a good reason to include a higher ebike share?

# Calculate the benefit of cycling and walking because they raise the MET minutes and thus lower
# the all cause moratality of low physical activity
# -> calculate the marginal benefit of one km cycling

#lowering of all cause mortality of one km 
# all cause mortality from here: http://ghdx.healthdata.org/record/ihme-data/gbd-2019-mortality-estimates-1950-2019

# low physical activity risk factor from here:
# http://ghdx.healthdata.org/gbd-results-tool?params=gbd-api-2019-permalink/f50245689196c894965621f48bf542fd


## linearize RR reduction of METs minutes per week

edge_internalization_mapping_low             = {"Cycle": "Cycling",
                                                "Walk": "Walking"}


# we need all cause mortality from all risk factors
df_low =  df.loc[(df['rei_name'] == 'All risk factors')& (df['measure_name'] == 'DALY') & (df['cause_name'] == 'All causes')]


# marginal benefit of all risk factors on all cause mortality according to \cite{Lear2017}.
# They give the effect of 0.85 RR for moderate (600-3000) vs low (<600) physical activity. 
# We use the middle of 1800 and get a marginal of 8.3*10-5 for one MET min per week.
# Which is 1.6*10-6 MET per year.

# calculate the benefit per MET
df_low['value'] = df_low['value'] * 1.6*10**-6





df_low['Variable'] = 'Cycling'
df_low['Unit'] = 'Physical activity benefit DALY/MET per capita'
df_low = df_low[['remind','year','Variable','Unit','value']]
df_low = df_low.rename(columns={"year":"Period", "remind": "Region",})
#add walking
df_low_walking = df_low.copy()
df_low_walking['Variable'] = 'Walking'


#multiply with MET_minutes for cycling, which is the only one available
df_low['value'] = df_low['value']*MET_minutes_cycling
df_low_walking['value'] = df_low_walking['value']*MET_minutes_walking
df_low = pd.concat([df_low,df_low_walking])