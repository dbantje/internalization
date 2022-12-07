from premise import NewDatabase
from premise.utils import eidb_label

from lca2rmnd.reporting import TransportLCAReporting
# from lca2rmnd.utils import project_string

import pandas as pd
import time
import glob
import re
import os
import gc
import brightway2 as bw
from pathlib import Path
import bw2io

# bw.projects.set_current("transport_lca")

remind_output_folder = "/home/alois/remind/testruns/lca_paper/v5/"
mifs = glob.glob(remind_output_folder + "*.mif")
storage_folder = "output/"
ldv_fleet_files = "/home/alois/git/premise/premise/data/iam_output_files/fleet_files/remind/passenger_cars/"

if "REMIND_generic_" in mifs[0]:
    for mif in mifs:
        os.rename(mif, mif.replace("REMIND_generic_", ""))

mifs = glob.glob(remind_output_folder + "*.mif")
scenarios = [
    re.search(
        '^remind_(.+?)\.mif',
        os.path.basename(txt)).group(1)
    for txt in mifs]

omitted_scens = [
]

scenarios = [sc for sc in scenarios if sc not in omitted_scens]

years = [
    2020, 2025, 2030, 2035, 2040, 2050
]

# scenarios = [
    # "Budg1100_H2Push_LowD",
    # "Budg1100_H2Push",
    # "Budg1100_ConvSyn_LowD",
    # "Budg1100_ConvSyn"
# ]
# years = [2050]

ecoinvent_version = "3.7.1"
eidb_original_label = 'ecoinvent {} cutoff'.format(ecoinvent_version)
fpei = "/home/alois/ecoinvent/ecoinvent {}_cutoff_ecoSpold02/datasets/".format(ecoinvent_version)

model = "remind"

# target recycling rates
# Pt: 30% to 80%, Co: 17-36%, Li-Ca: 13-28%
recycling_paradigm = "high"
target_rates={
    "high": {
        "platinum": 0.8, "cobalt": 0.36, "lithium carbonate": 0.28
    },
    "low": {
        "platinum": 0.3, "cobalt": 0.17, "lithium carbonate": 0.13
    }
}

scencars = [
    "transport, passenger car, fleet average, fuel cell electric",
    "transport, passenger car, fleet average, gasoline",
    "transport, passenger car, fleet average, battery electric",
    "transport, passenger car, fleet average, gasoline",
]

mp2eps = {
    "human health": [
        ('ReCiPe Midpoint (H) V1.13', 'human toxicity', 'HTPinf'),
        # ('ReCiPe Midpoint (H) V1.13', 'climate change', 'GWP100'),
        ('ReCiPe Midpoint (H) V1.13', 'ionising radiation', 'IRP_HE'),
        ('ReCiPe Midpoint (H) V1.13', 'ozone depletion', 'ODPinf'),
        ('ReCiPe Midpoint (H) V1.13', 'particulate matter formation', 'PMFP'),
        ('ReCiPe Midpoint (H) V1.13', 'photochemical oxidant formation', 'POFP'),
    ],
    "ecosystem quality": [
        ('ReCiPe Midpoint (H) V1.13', 'freshwater ecotoxicity', 'FETPinf'),
        ('ReCiPe Midpoint (H) V1.13', 'marine ecotoxicity', 'METPinf'),
        ('ReCiPe Midpoint (H) V1.13', 'terrestrial ecotoxicity', 'TETPinf'),
        ('ReCiPe Midpoint (H) V1.13', 'agricultural land occupation', 'ALOP'),
        # ('ReCiPe Midpoint (H) V1.13', 'climate change', 'GWP100'),
        ('ReCiPe Midpoint (H) V1.13', 'freshwater eutrophication', 'FEP'),
        ('ReCiPe Midpoint (H) V1.13', 'marine eutrophication', 'MEP'),
        # ('ReCiPe Midpoint (H) V1.13', 'natural land transformation', 'NLTP'),
        ('ReCiPe Midpoint (H) V1.13', 'terrestrial acidification', 'TAP100'),
        ('ReCiPe Midpoint (H) V1.13', 'urban land occupation', 'ULOP'),
    ],
    "resources": [
        ('ReCiPe Midpoint (H) V1.13', 'metal depletion', 'MDP'),
        ('ReCiPe Midpoint (H) V1.13', 'fossil depletion', 'FDP'),
        ('ReCiPe Midpoint (H) V1.13', 'water depletion', 'WDP')
    ]
}

co2met = ('Custom Method', 'CO2 emissions', 'kgCO2')


def project_string(scenario, project="transport"):
    return "{}_lca_{}".format(project, scenario)

def create_endpoint_methods():
    """
    Create mid2endpoint methods starting from all midpoint methods
    and extracting relevant CFs from the endpoint total methods.
    """

    if len([m for m in bw.methods if m[0] == "ReCiPe Endpoint (H,A) V1.13"]) == 0:
        bw2io.package.BW2Package.import_file("/home/alois/git/lca_paper/ReCiPe-2008-v113-endpoint.bw2package")

    units = {
        "human health": "DALY",
        "ecosystem quality": "species lost.year",
        "resources": "USD2005"
    }

    ep_categories = {
        "human health": "Human health",
        "ecosystem quality": "Ecosystems",
        "resources": "Resources"
    }

    # remove existing ep methods
    eps = [m for m in bw.methods if m[0] == 'ReCiPe Endpoint (H) V1.13']
    for m in eps:
        bw.Method(m).deregister()

    for cat, ep in ep_categories.items():
        method = ('ReCiPe Endpoint (H,A) V1.13', ep)
        odm = bw.Method(method)
        fcts = odm.load()
        eps = pd.Series({(f[0], bw.get_activity(f[0])["name"]): f[1] for f in fcts})\
                .rename_axis(["key", "flow"]).to_frame()
        midpoint_methods = [m for m in bw.methods if m[0] == 'ReCiPe Midpoint (H) V1.13']

        for mm in mp2eps[cat]:
            odm = bw.Method(mm)
            fcts = odm.load()
            mps = pd.Series({(f[0], bw.get_activity(f[0])["name"]): f[1] for f in fcts})\
                   .rename_axis(["key", "flow"]).to_frame()
            cfs = mps.merge(eps, left_index=True, right_index=True)
            if len(cfs) == 0:
                print("No CFs found in Endpoint total for {}".format(mm))
                continue
            dct = cfs.reset_index("flow")["0_y"].to_dict()

            newfcts = [(key, value) for key, value in dct.items()]
            # save to new EP method
            mname = ('ReCiPe Endpoint (H) V1.13', cat, mm[1])
            epm = odm.copy(mname)
            epm.register(unit=units[cat])
            epm.write(newfcts)


def fix_water_usage_method():
    for scenario in scenarios:
        bw.projects.set_current(project_string(scenario))
        print("Fixing water usage method for {}".format(scenario))
        method = ('ReCiPe Midpoint (H) (obsolete)', 'water depletion', 'WDP')
        odm = bw.Method(method)
        fcts = odm.load()
        newfcts = [(el[0], 1) for el in fcts]
        odm.write(newfcts)

def add_endpoint_methods():
    for scenario in scenarios:
        bw.projects.set_current(project_string(scenario))
        print("Add enpoint methods to {}".format(scenario))
        create_endpoint_methods()

def create_co2_method():
    print("Creating CO2 methods for {}.".format(bw.projects.current))
    cm = bw.Method(('IPCC 2013', 'climate change', 'GWP 100a'))

    fcts = cm.load()
    mps = pd.Series({(f[0], bw.get_activity(f[0])["name"]): f[1] for f in fcts})\
            .rename_axis(["key", "flow"]).to_frame().reset_index("flow")
    mps = mps.loc[mps.flow.str.startswith("Carbon dioxide")]
    dct = mps[0].to_dict()
    newfcts = [(key, value) for key, value in dct.items()]

    com = cm.copy(co2met)
    com.register(unit="kg CO2")
    com.write(newfcts)

def add_co2method():
    for scenario in scenarios:
        bw.projects.set_current(project_string(scenario))
        if len([m for m in bw.methods if m == co2met]) == 0:
            create_co2_method()

def add_recipe2016():
    for scenario in scenarios:
        bw.projects.set_current(project_string(scenario))
        print("Add new methods to {}".format(project_string(scenario)))
        # remove existing ep methods
        eps = [m for m in bw.methods if m == ('ReCiPe Midpoint (H) V1.13', 'water consumption', "WCO")]
        for m in eps:
            bw.Method(m).deregister()

        bw2io.package.BW2Package.import_file("/home/alois/git/lca_paper/ReCiPe-2016.bw2package")

        wm = bw.Method(('ReCiPe 2016', '1.1 (20180117)', 'Midpoint', 'Water consumption'))
        new_wm = wm.copy(('ReCiPe Midpoint (H) V1.13', 'water consumption', "WCO"))
        new_wm.register()
        wm.deregister()

def scenario_setup():
    if "biosphere3" not in bw.databases:
        print("Setting up missing biosphere.")
        bw.bw2setup()
    if eidb_original_label not in bw.databases:
        print("Importing ecoinvent DB.")
        ei = bw.SingleOutputEcospold2Importer(
            fpei, eidb_original_label)
        ei.apply_strategies()
        ei.statistics()
        ei.write_database()
    # create_endpoint_methods()
    # create_co2_method()
    # fix_water_usage_method()

def create_inventories():
    """
    Create REMIND-ecoinvent inventories for given years and
    scenarios. Inventories for scenarios that are to be
    omitted are deleted from the database to speed up
    the database.
    """

    for scenario in scenarios:
        print("Creating inventories for {}.".format(scenario))
        pstr = project_string(scenario)
        bw.projects.set_current(pstr)
        scenlist = []
        dblist = [db for db in bw.databases if db.startswith("ecoinvent_{}".format(model))]

        for year in years:
            if "ecoinvent_{}_{}_{}".format(model, scenario, year) in dblist:
                print("Existing database found for {}, year {}.".format(scenario, year))
                continue
            scen = {"model": model, "pathway": scenario, "year": year,
                    "passenger_cars": {
                        "regions":["EUR"],
                        "fleet file": os.path.join(
                            ldv_fleet_files,
                            "_".join([scenario, "vintcomp.csv"]))},
                    "filepath": remind_output_folder,
            }
            if year <= 2024:
                scen["exclude"] = ["update_electricity"]
            # scenlist.append(scen)
            scenario_setup()
            start_time = time.time()
            ndb = NewDatabase(scenarios=[scen],
                              source_db=eidb_original_label,
                              source_version=ecoinvent_version)

            ndb.update_all()
            ndb.write_db_to_brightway()
            print("Creation of scenario {} took: {:.1f} min"
                  .format(scenario, (time.time() - start_time)/60.))
            del ndb
            gc.collect()

        # if not scenlist:
        #     print("All inventories present for {}".format(scenario))
        #     continue


def purge():
    for scenario in scenarios:
        if project_string(scenario) in [p.name for p in bw.projects]:
            print("Purging project {}".format(project_string(scenario)))
            bw.projects.delete_project(project_string(scenario), delete_dir=True)


def clean():
    for scenario in scenarios:
        bw.projects.set_current(project_string(scenario))
        print("Cleaning scenario {}".format(scenario))
        invlist = list(bw.databases)
        for db in invlist:
            if db.startswith("ecoinvent_"):
                print("Deleting {}".format(db))
                del bw.databases[db]


def get_storage_file_path(scenario, category, rec_paradigm=None):
    if rec_paradigm is None:
        rec_paradigm = recycling_paradigm
    rec_subfolder = "recycling_{}".format(rec_paradigm)
    subpath = os.path.join(storage_folder, rec_subfolder)
    if not os.path.exists(subpath):
        os.mkdir(subpath)
    return os.path.join(subpath, "{}_{}.csv".format(category, scenario))

def fix_electricity_markets():
    region = "EUR"
    for scenario in scenarios:
        print("Fixing electricity markets for {}".format(scenario))
        bw.projects.set_current(project_string(scenario))
        for year in years:
            db = bw.Database(eidb_label(model, scenario, year))
            # if year > 2020:
            elec_market = [a for a in db if a["name"] == "market group for electricity, low voltage" and a["location"] == "EUR"][0]
            # else:
            #     continue
            # fix for BEVs and PHEVs
            fleets = [
                "transport, passenger car, fleet average, plugin gasoline hybrid",
                "transport, passenger car, fleet average, plugin diesel hybrid",
                "transport, passenger car, fleet average, battery electric"
            ]
            for fname in fleets:
                fleet_act = [a for a in db if a["name"] == fname and a["location"] == region][0]
                exc = [ex for ex in fleet_act.technosphere() if ex["product"] == "electricity, low voltage"]
                assert len(exc) == 1
                exc = exc[0]
                exc.input = elec_market
                exc.save()

def fix_bev_electricity_markets():
    for scenario in scenarios:
        print("Fixing BEV electricity markets for {}".format(scenario))
        bw.projects.set_current(project_string(scenario))
        db = bw.Database(eidb_label(model, scenario, 2020))
        supply = [a for a in db if a["name"] == "electricity supply for electric vehicles" and a["location"] == "EUR"][0]
        market = [a for a in db if a["name"] == "market group for electricity, low voltage" and a["location"] == "RER"][0]

        # replace existing
        exc = [ex for ex in supply.technosphere()]
        assert len(exc) == 1
        exc = exc[0]
        exc.input = market
        exc.save()

# for scenario in scenarios:
#     if "H2Push" in scenario:
#         print("Fixing hydrogen markets for {}".format(scenario))
#         bw.projects.set_current(project_string(scenario))
#         db = bw.Database(eidb_label(model, scenario, 2020))
#         smr_name = "Hydrogen, gaseous, 700 bar, from SMR of NG, at fuelling station"
#         # replace supply names in case!
#         elec_name = "Hydrogen, gaseous, 700 bar, from electrolysis, at fuelling station"

#         smr_act = [a for a in db if a["name"] == smr_name and a["location"] == "RER"][0]
#         supply_act = [a for a in db if a["location"] == region and a["name"] == "fuel supply for hydrogen vehicles"][0]

#         # old exchange:
#         elec_ex = [ex for ex in supply_act.technosphere()]
#         for ex in elec_ex:
#             ex.delete()
#         smr_ex = supply_act.new_exchange(input=smr_act, amount=1, type="technosphere")
#         smr_ex.save()

def fix_hydrogen_markets():
    region = "EUR"
    # load shares
    data = pd.concat((
        pd.read_csv(
            os.path.join(remind_output_folder, "remind_{}.mif".format(sc)),
            sep=";", index_col=["Region", "Variable"])
        .loc[region]
        .assign(Scenario=sc)
        for sc in scenarios)).drop(columns=["Model", "Unit", "Unnamed: 24"])
    variables = [
        "SE|Hydrogen|Electricity",
        "SE|Hydrogen"
    ]

    data = data.loc[variables].reset_index()
    data = data.melt(id_vars=["Variable", "Scenario"], var_name="Year")\
               .pivot(index=["Scenario", "Year"], columns="Variable", values="value")
    data["share"] = data["SE|Hydrogen|Electricity"] / data["SE|Hydrogen"]

    for scenario in scenarios:
        if "ConvSyn" in scenario:
            print("Fixing hydrogen markets for {}".format(scenario))
            bw.projects.set_current(project_string(scenario))
            for year in years:
                db = bw.Database(eidb_label(model, scenario, year))
                smr_names = {
                    "gasoline": "Gasoline, synthetic, from MTG, hydrogen from SMR of natural gas, energy allocation, at fuelling station",
                    "diesel": "Diesel, synthetic, from natural gas-based hydrogen, energy allocation, at fuelling station"
                }
                # replace supply names in case!
                elec_names = {
                    "diesel": "Diesel, synthetic, from electrolysis-based hydrogen, energy allocation, at fuelling station",
                    "gasoline": "Gasoline, synthetic, from MTG, hydrogen from electrolysis, energy allocation, at fuelling station"
                }
                for fuel, smr_supply in smr_names.items():
                    smr_act = [a for a in db if a["name"] == smr_supply and a["location"] == "RER"][0]
                    supply_act = [a for a in db if a["location"] == region and a["name"] == "fuel supply for {} vehicles".format(fuel)][0]
                    elec_act = [a for a in db if a["name"] == elec_names[fuel] and a["location"] == "RER"][0]

                    # old exchange:
                    elec_ex = [ex for ex in supply_act.technosphere() if ex.input == elec_act]
                    if len(elec_ex) == 0:
                        print("No {} supply found for {}, {}".format(elec_act["name"], scenario, year))
                        continue
                    elec_ex = elec_ex[0]
                    smr_amount = elec_ex["amount"] * (1-data.loc[(scenario, str(year)), "share"])
                    smr_ex = supply_act.new_exchange(input=smr_act, amount=smr_amount, type="technosphere")

                    elec_ex["amount"] *= data.loc[(scenario, str(year)), "share"]
                    smr_ex.save()
                    elec_ex.save()
        if "H2Push" in scenario:
            print("Fixing hydrogen markets for {}".format(scenario))
            bw.projects.set_current(project_string(scenario))
            for year in years:
                db = bw.Database(eidb_label(model, scenario, year))
                smr_name = "Hydrogen, gaseous, 700 bar, from SMR of NG, at fuelling station"
                # replace supply names in case!
                elec_name = "Hydrogen, gaseous, 700 bar, from electrolysis, at fuelling station"

                smr_act = [a for a in db if a["name"] == smr_name and a["location"] == "RER"][0]
                supply_act = [a for a in db if a["location"] == region and a["name"] == "fuel supply for hydrogen vehicles"][0]
                elec_act = [a for a in db if a["name"] == elec_name and a["location"] == "RER"][0]

                # old exchange:
                elec_ex = [ex for ex in supply_act.technosphere() if ex.input == elec_act]
                if len(elec_ex) == 0:
                    print("No {} supply found for {}, {}".format(elec_act["name"], scenario, year))
                    continue
                elec_ex = elec_ex[0]
                if year == 2020:
                    elec_ex.input = smr_act
                    elec_ex.save()
                else:
                    smr_amount = elec_ex["amount"] * (1-data.loc[(scenario, str(year)), "share"])
                    smr_ex = supply_act.new_exchange(input=smr_act, amount=smr_amount, type="technosphere")

                    elec_ex["amount"] *= data.loc[(scenario, str(year)), "share"]
                    smr_ex.save()
                    elec_ex.save()



def fix_tags():
    for scenario in scenarios:
        print("Fixing tags for {}".format(scenario))
        bw.projects.set_current(project_string(scenario))
        for year in years:
            db = bw.Database(eidb_label(model, scenario, year))
            for car in scencars:
                act = [a for a in db if a["name"] == car][0]
                for ex in act.exchanges():
                    if "tag" not in ex:
                        print("No tag in {}, {}".format(ex, db.name))
                        if ex["name"].startswith("market group for electricity"):
                            ex["tag"] = "energy chain"
                            ex.save()


def restore():
    for scenario in scenarios:
        backup = glob.glob("ei_data/brightway2-project-{}-backup.*.tar.gz".format(project_string(scenario)))
        assert len(backup) == 1
        print("Restoring {}".format(scenario))
        bw.restore_project_directory(backup[0])
    update_recycling_rates()
    fix_tags()
    add_co2method()
    fix_hydrogen_markets()
    fix_electricity_markets()


def inventory_check():
    projectlst = [p.name for p in bw.projects
                  if p.name != "default"
                  and not p.name.startswith("REMIND_")]
    for pro in projectlst:
        print("Project: {}".format(pro))
        bw.projects.set_current(pro)
        for db in bw.databases:
            if db not in ["biosphere3", eidb_original_label]:
                print("{}: {}".format(db, len(bw.Database(db))))


def update_recycling_rates():
    print("Updating recycling rates, paradigm: {}".format(recycling_paradigm))
    bw.projects.set_current(project_string(scenarios[0]))
    # obtain original rates
    orig_db = bw.Database(eidb_original_label)
    orig_rates = {}
    for material in target_rates["high"]:
        act = [a for a in orig_db if a["name"] == "market for {}".format(material)]
        assert len(act) == 1
        act = act[0]
        sum_recyc = 0
        for ex in act.technosphere():
            if material == ex["name"] and "treatment" in bw.get_activity(ex["input"])["name"]:
                sum_recyc += ex["amount"]
        sum_prod = 1-sum_recyc
        orig_rates[material] = sum_recyc/sum_prod

    for scenario in scenarios:
        bw.projects.set_current(project_string(scenario))
        print("Updating rates for {}".format(scenario))
        for year in years:
            db = bw.Database(eidb_label(model, scenario, year))
            rec_rates = target_rates[recycling_paradigm]
            for material in rec_rates:
                act = [a for a in db if a["name"] == "market for {}".format(material)]
                assert len(act) == 1, "None or too many activities found for {} in {}".format(material, db.name)
                act = act[0]
                # primary production
                prod_exs = [ex for ex in act.technosphere()
                            if "treatment" not in ex["name"]
                            and material == ex["product"]]
                tot_prod = sum([ex["amount"] for ex in prod_exs])
                recyc_exs = [ex for ex in act.technosphere() if "treatment" in ex["name"]
                             and material == ex["name"]]
                tot_recyc = sum([ex["amount"] for ex in recyc_exs])

                assert year >= 2020
                new_rate = orig_rates[material] + (rec_rates[material]-tot_recyc) * min(1, (year - 2020)/30)
                for ex in prod_exs:
                    ex["amount"] = (1-new_rate) * ex["amount"]/tot_prod
                    ex.save()
                for ex in recyc_exs:
                    ex["amount"] = new_rate * ex["amount"]/tot_recyc
                    ex.save()

def restore_recycling_rates():
    for scenario in scenarios:
        print("Restoring rates for {}".format(scenario))
        bw.projects.set_current(project_string(scenario))
        db = bw.Database(eidb_original_label)
        original_technosphere = {}
        for material in target_rates["high"]:
            act = [a for a in db if a["name"] == "market for {}".format(material)]
            assert len(act) == 1
            # primary production
            original_technosphere[material] = {
                bw.get_activity(ex["input"])["name"]: ex["amount"] for ex in act[0].technosphere()
                if ex["type"] == "technosphere" and material == ex["name"]}

        for year in years:
            db = bw.Database(eidb_label(model, scenario, year))
            for material in materials:
                act = [a for a in db if a["name"] == "market for {}".format(material)]
                assert len(act) == 1
                for ex in act[0].technosphere():
                    if ex["type"] == "technosphere" and material == ex["product"]:
                        ex["amount"] = original_technosphere[material][bw.get_activity(ex["input"])["name"]]
                        ex.save()

def export():
    for scenario in scenarios:
        bw.backup_project_directory(project_string(scenario))

def store_output():
    bw.projects.set_current(project_string(scenarios[0]))

    if not os.path.exists(storage_folder):
        os.mkdir(storage_folder)

    for scenario in scenarios:
        start_time = time.time()
        print("Reporting LCA data for {}".format(scenario))
        omit = ["climate change", "natural land transformation"]
        indicators = [m for m in bw.methods if m[0] == "ReCiPe Endpoint (H) V1.13"
                      and m[2] not in omit]
        indicators.extend([
            ("IPCC 2013", "climate change", "GWP 100a"),
            co2met])
        rep = TransportLCAReporting(
            scenario, years, project_string(scenario), Path(remind_output_folder),
            methods=indicators,
            regions=["EUR"])
        # print("Skipping calculation of LCA results!")
        ldv_lca = rep.report_LDV_LCA()
        ldv_lca.to_csv(get_storage_file_path(scenario, "LDV_LCA"))
        materials = rep.report_materials()
        materials.to_csv(get_storage_file_path(scenario, "materials"))
        # direct_emis = rep.report_direct_emissions()
        # direct_emis.to_csv(get_storage_file_path(scenario, "direct_emissions"))
        print("Report for scenario {} took: {:.1f} min"
              .format(scenario, (time.time() - start_time)/60.))

# fix BEV electricity supply
# for scenario in scenarios:
#     bw.projects.set_current(project_string(scenario))
#     print("Updating electricity supply for {}".format(scenario))
#     for year in years:
#         db = bw.Database(eidb_label(model, scenario, year))
#         bev = [a for a in db if a["name"] == "transport, passenger car, fleet average, battery electric"][0]
#         supply = [a for a in db if a["name"] == "market group for electricity, low voltage"
#                   and a["location"] == "EUR"][0]
#         bev_ex = [
#                 e for e in bev.technosphere() if e["name"] == "market group for electricity, low voltage"
#         ]
#         assert len(bev_ex) == 1
#         bev_ex = bev_ex[0]
#         bev_ex.input = supply
#         bev_ex.save()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            'Prepare hybrid REMIND/ecoinvent inventories, '
            'store output and produce overview plots.'))
    parser.add_argument('--remind_folder', metavar='path', required=False,
                        help='the path to the REMIND output folder')
    parser.add_argument('--storage_folder', metavar='path', required=False,
                        help='path for output data storage')
    parser.add_argument('--plot-folder', metavar='path', required=False,
                        help='path for output data storage')
    parser.add_argument('-c', dest='create', action='store_true',
                        help='create inventories')
    parser.add_argument('-s', dest='store', action='store_true',
                        help='store reporting data')
    parser.add_argument('--clean', action='store_true',
                        help='remove existing inventories')
    parser.add_argument('--recreate', action='store_true',
                        help='clean and recreate inventories')
    parser.add_argument('--update-recycling', action='store_true',
                        help='update recycling rates for batteries and fuel cells')
    parser.add_argument('--restore-recycling', action='store_true',
                        help='restore recycling rates for batteries and fuel cells')
    parser.add_argument('--fix-tags', action='store_true',
                        help='fix missing energy chain tags')
    parser.add_argument('--fix-hydrogen-markets', action='store_true',
                        help='fix hydrogen markets')
    parser.add_argument('--fix-electricity-markets', action='store_true',
                        help='fix electricity markets')
    parser.add_argument('--fix-bev-electricity-markets', action='store_true',
                        help='fix BEV electricity markets')
    parser.add_argument('--add-co2', action='store_true',
                        help='add CO2 method to scenarios')
    parser.add_argument('--add-recipe', action='store_true',
                        help='add ReCiPe2016 methods to scenarios')
    parser.add_argument('--add-endpoints', action='store_true',
                        help='add updated endpoint methods to scenarios')
    parser.add_argument('--fix-water', action='store_true',
                        help='fix water usage method.')
    parser.add_argument('--purge', action='store_true',
                        help='remove all inventories (including biosphere and ecoinvent)')
    parser.add_argument('--restore', action='store_true',
                        help='restore project directories')
    parser.add_argument('--full', action='store_true',
                        help='clean and recreate inventories, generate impact tables and plot results')
    parser.add_argument('--setup', action='store_true',
                        help='create scenario specific projects and import biosphere and technosphere')
    parser.add_argument('--export', action='store_true',
                        help='export existing scenario project directories')
    args = parser.parse_args()

    if args.remind_folder is not None:
        remind_output_folder = args.remind_folder
    if args.plot_folder is not None:
        plot_folder = args.plot_folder
    if args.storage_folder is not None:
        storage_folder = args.storage_folder
    if args.create:
        create_inventories()
    if args.store:
        store_output()
    if args.clean:
        clean()
    if args.update_recycling:
        update_recycling_rates()
    if args.restore_recycling:
        restore_recycling_rates()
    if args.fix_tags:
        fix_tags()
    if args.fix_hydrogen_markets:
        fix_hydrogen_markets()
    if args.fix_electricity_markets:
        fix_electricity_markets()
    if args.fix_bev_electricity_markets:
        fix_bev_electricity_markets()
    if args.add_co2:
        add_co2method()
    if args.add_recipe:
        add_recipe2016()
    if args.add_endpoints:
        add_endpoint_methods()
    if args.fix_water:
        fix_water_usage_method()
    if args.purge:
        purge()
    if args.restore:
        restore()
    if args.recreate:
        clean()
        create_inventories()
    if args.full:
        clean()
        create_inventories()
        store_output()
    if args.setup:
        scenario_setup()
    if args.export:
        export()
