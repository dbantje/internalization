
# load the GBD data and join it with mappings
with open(GBD_path,'r') as f:
    data = json.loads(f.read())
    
# Normalizing data
df = pd.json_normalize(data, record_path =['data'])
cols = pd.json_normalize(data, record_path =['cols'])
df.columns = cols[0]


# sex is boths, age is all - can be deleted
df = df.drop(['sex', 'age'], axis=1)

# join to the mappings for measure, location, cause, metric and rei
gbd_location_mapping = pd.read_csv(gbd_location_mapping_path, sep = ';', encoding='latin-1')
gbd_measure_mapping  = pd.read_csv(gbd_measure_mapping_path , encoding='latin-1')
gbd_rei_mapping      = pd.read_csv(gbd_rei_mapping_path, sep = ';', encoding='latin-1' )
gbd_cause_mapping    = pd.read_csv(gbd_cause_mapping_path, sep = ';'   , encoding='latin-1')
gbd_metric_mapping   = pd.read_csv(gbd_metric_mapping_path  , encoding='latin-1')

df = df.merge(gbd_location_mapping,  how='left')
df = df.merge(gbd_measure_mapping,  how='left')
df = df.merge(gbd_rei_mapping,  how='left')
df = df.merge(gbd_cause_mapping,  how='left')
df = df.merge(gbd_metric_mapping,  how='left')
df = df.drop(['measure', 'location','cause','rei','metric','upper','lower'], axis=1)

df = df.pivot_table(index=['year','location_name','remind','measure_name','rei_name','cause_name'], 
                    columns='metric_name', 
                    values='val')


# aggregate rate to remind regions with weighted average of cases
df = df.groupby(["year", "remind","measure_name","rei_name","cause_name"]).apply(wavg, "rate", "number")
df = pd.DataFrame(df)

df = df.reset_index()
df['value'] = df[0] / 100000
df = df.drop([0], axis=1)

