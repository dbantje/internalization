#read nutrition intake and supply from MAGPIE results
#calculate the shares of the total intake and total supply
magpie = pd.read_csv(magpie_report, sep=';')
magpie = magpie[['Model', 'Scenario', 'Region', 'Variable', 'Unit',str(year)]]
magpie_intake = magpie.query('Variable.str.contains("Nutrition|Calorie Intake") and Unit.str.contains("kcal/capita/day") and Region.str.contains("World")', engine='python')
magpie_intake_total = magpie.query('Variable == ("Nutrition|Calorie Intake") and Unit.str.contains("kcal/capita/day") and Region.str.contains("World")', engine='python')
magpie_supply = magpie.query('Variable.str.contains("Nutrition|Calorie Supply") and Unit.str.contains("kcal/capita/day") and Region.str.contains("World")', engine='python')
magpie_supply_total = magpie.query('Variable == ("Nutrition|Calorie Supply") and Unit.str.contains("kcal/capita/day") and Region.str.contains("World")', engine='python')

#shares
magpie_share_supply_intake = magpie_supply[str(year)].div(magpie_intake[str(year)])
magpie_intake[str(year)] = magpie_intake[str(year)].div(int(magpie_intake_total[str(year)]))

#intake scaled with supply to account for waste
#this introduces a small error leading to higher impacts since meat has a lower food waste than the total average
magpie_intake[str(year)] = magpie_intake[str(year)]* magpie_share_supply_intake
magpie_intake['Unit'] = "share of kcal nutrition supply, normalized by demand"


#select which nutrition processes are included
#map ecoinvent processes to the MAGPIE variables
# tropical root and eggs are mapped to the closest proces
magpie_ecoinvent_nutrition_mapping_extended = {'Nutrition|Calorie Intake|+|Fish': "'market for fishmeal, 65-67% protein' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Crops|Cereals|+|Maize': "'market for maize seed, organic, at farm' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Crops|Cereals|+|Rice': "'market for rice, non-basmati' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Crops|Cereals|+|Temperate cereals': "'market for wheat grain, feed' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Crops|+|Oil crops': "'market for peanut' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Crops|Other crops|+|Fruits Vegetables Nuts': "'market for apple' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Crops|Other crops|+|Potatoes': "'market for potato' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Crops|Other crops|+|Pulses': "'market for chickpea' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Crops|Other crops|+|Tropical roots': "'market for potato' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Crops|Sugar crops|+|Sugar beet': "'market for sugar beet pulp' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Crops|Sugar crops|+|Sugar cane': "'market for sugar, from sugarcane' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Livestock products|+|Dairy': "'market for dairy' (cubic meter, GLO, None)",
                                               'Nutrition|Calorie Intake|Livestock products|+|Eggs': "'chicken production' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Livestock products|+|Monogastric meat': "'market for swine for slaughtering, live weight' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Livestock products|+|Poultry meat': "'chicken production' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Livestock products|+|Ruminant meat': "'market for red meat, live weight' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Secondary products|+|Alcoholic beverages': "'market for ethoxylated alcohol (AE11)' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Secondary products|+|Brans': "'market for wheat grain, feed' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Secondary products|+|Molasses': "'market for molasses, from sugar beet' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Secondary products|+|Oils': "'market for palm oil, refined' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Secondary products|+|Sugar': "'market for sugar, from sugarcane' (kilogram, GLO, None)",
                                               'Nutrition|Calorie Intake|Crops|Cereals|+|Tropical cereals': "'market for millet' (kilogram, GLO, None)"}
											   
											   
#kg to kcal
#source USDA National Nutrient Database
ecoinvent_kcal_kg_mapping_extended          = {"'market for fishmeal, 65-67% protein' (kilogram, GLO, None)": 3000,
                                               "'market for maize seed, organic, at farm' (kilogram, GLO, None)": 3840,
                                               "'market for rice, non-basmati' (kilogram, GLO, None)": 1300,
                                               "'market for wheat grain, feed' (kilogram, GLO, None)": 3400,
                                               "'market for peanut' (kilogram, GLO, None)": 5670,
                                               "'market for apple' (kilogram, GLO, None)": 520,
                                               "'market for potato' (kilogram, GLO, None)": 770,
                                               "'market for chickpea' (kilogram, GLO, None)": 3780,
                                               "'market for potato' (kilogram, GLO, None)": 770,
                                               "'market for sugar beet pulp' (kilogram, GLO, None)": 430,
                                               "'market for sugar, from sugarcane' (kilogram, GLO, None)": 760,
                                               "'market for dairy' (cubic meter, GLO, None)": 460000,
                                               "'chicken production' (kilogram, GLO, None)": 1430,
                                               "'market for swine for slaughtering, live weight' (kilogram, GLO, None)": 2600,
                                               "'chicken production' (kilogram, GLO, None)": 1430,
                                               "'market for red meat, live weight' (kilogram, GLO, None)": 1980,
                                               "'market for ethoxylated alcohol (AE11)' (kilogram, GLO, None)": 2500,
                                               "'market for wheat grain, feed' (kilogram, GLO, None)": 3400,
                                               "'market for molasses, from sugar beet' (kilogram, GLO, None)": 2900,
                                               "'market for palm oil, refined' (kilogram, GLO, None)": 8840,
                                               "'market for sugar, from sugarcane' (kilogram, GLO, None)": 3800,
                                               "'market for millet' (kilogram, GLO, None)": 3780,}
											   
magpie_intake['Variable'] = magpie_intake['Variable'].map(magpie_ecoinvent_nutrition_mapping_extended)

#calculate the kg/kcal required food intake of the different foods
magpie_intake['kcal/kg'] = magpie_intake['Variable']
magpie_intake['kcal/kg'] = magpie_intake['kcal/kg'].map(ecoinvent_kcal_kg_mapping_extended)
magpie_intake['kcal/kg'] = 1 / magpie_intake['kcal/kg']
magpie_intake[str(year)] = magpie_intake[str(year)] * magpie_intake['kcal/kg']
magpie_intake['Unit'] = "share of kcal nutrition supply per one kcal, normalized by demand"
magpie_intake = magpie_intake.drop(['kcal/kg'], axis=1)


