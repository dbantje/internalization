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

# In this notebook, we create a new project, and import ecoinvent 3.5 with the cutoff system model.

# Create new project and import libraries:

#projects

#projects.current

# Check current version of important python packages. This code was written with Brightway 2.1.1, bw2io 0.6 RC2 and wurst 0.1.
#
# update: Oct 31, 2018:
#
# brightway2 (2, 3)
#
# bw2io (0, 7, 1)
#
# wurst (0, 1, 2)
#

# +
from brightway2 import *
#print ('brightway2', brightway2.__version__)


from bw2data.project import ProjectDataset	


import bw2io
print ('bw2io', bw2io.__version__)

import wurst
print ('wurst', wurst.__version__)

# -

# Prepare brightway:

bw2setup()

# Import ecoinvent:

bw2io.add_ecoinvent_35_biosphere_flows()

import os
from os.path import *
ecoinvent_filepath =  os.path.expanduser("data/ecoinvent/datasets")
brightway_filepath = os.path.expanduser("data/brightway/")


#switch directory path

assert os.path.isdir(brightway_filepath)
projects._base_data_dir = brightway_filepath
projects._base_logs_dir = brightway_filepath+ "logs"
#os.mkdir(projects._base_logs_dir)


ei = SingleOutputEcospold2Importer(ecoinvent_filepath, "ecoinvent_3.5")
ei.apply_strategies()
ei.statistics()
ei.write_database()

# +
#ei.drop_unlinked(i_am_reckless=True)
# -

databases

# %run initialize_notebook.ipynb
from helpers import lcahelp, eimod, ei2rmnd

# +
# Prepare additional datasets

#We import some additional datasets for carbon capture and storage electricity. These have been prepared in Simapro and exported. We have cleaned up the datasets and stored them in excel using Brightway input output functionality. We read them from excel here:

## Import Carma CCS database from excel


lcahelp.import_karma()


# +
# Define REMIND scenarios and years to be used:

#Here we set the scenarios and years for which we want to use when creating a new version of ecoinvent:

database_dict = {}

for year in years:
#for year in [2015]:
    for scenario in scenarios:
        if year == 2015 and scenario != 'BAU':
            continue
        db_name = 'ecoinvent_Remind_' + scenario + '_' + str(year)
        database_dict[db_name] = {'year' : year , 'scenario' : scenario}

#database_dict['test'] = {'year' : 2041 , 'scenario' : 'BAU'}
#for scenario in ['Base','NDC','NPi','PkBudg900','PkBudg1100','PkBudg1300']:

# +

database_dict


# +
databases



# +
# Modify all electricity generation and market datasets:

#This cell uses the Wurst functionality as well as the functions that we wrote in the notebook "_Functions_to_modify_ecoinvent.ipynb" to make the changes to ecoinvent:

eimod.modify_electricity_generation_datasets(database_dict)



# -

import brightway2 as bw

bw.projects.set_current('default')

#change so it write to the brightway folder
# save name and give it to next script

backup_project_directory('default')

#save name to use for next scroüpt

bw.projects

