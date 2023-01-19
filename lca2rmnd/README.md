# lca2rmnd
Calculate LCA impacts from REMIND output files

## Usage

To generate a reporting for the REMIND power sector, use the class `ElectricityLCAReporting` like so:
```
        rep = ElectricityLCAReporting(
            scenario, years, bw_project, Path(remind_output_folder),
            methods=indicators,
            regions=["EUR"])
        el_lca = rep.report_sectoral_LCA()
        el_lca.to_csv("sectoral_output.csv")
        
        tech_lca = rep.report_tech_LCA()

```
for the REMIND scenario `scenario` and the timesteps `years`, `bw_project` is the brightway2 project that hosts the relevant databases, 
`remind_output_folder` contains the REMIND scenario outputs, `methods` is a list of brightway2 methods for which one wants to generate reporting variables
and with `regions` the regional scope can be defined (has to be a valid REMIND region and part of the outputs, of course).

`report_sectoral_LCA` returns pandas tables with the year, region, (lca) method dimensions. Results are given per kWh and total.

 `report_tech_LCA` tries to find a set of activities for a REMIND output variable in the region and calculates impacts on the level of REMIND
 output variables. This tool has not been used for some time.
 
 
 Reporting for the transport sector can be done in a similar way:
 ```
        rep = TransportLCAReporting(
            scenario, years, project_string(scenario), Path(remind_output_folder),
            methods=indicators,
            regions=["EUR"])
        ldv_lca = rep.report_LDV_LCA()
        ldv_lca.to_csv(get_storage_file_path(scenario, "LDV_LCA"))
        materials = rep.report_materials()
        materials.to_csv(get_storage_file_path(scenario, "materials"))
```
`report_LDV_LCA` reports impacts for the REMIND energy service output variables for LDVs.
`report_materials` sums up the amounts of raw materials used (biosphere exchanges) for a set of materials that is available in ecoinvent.
